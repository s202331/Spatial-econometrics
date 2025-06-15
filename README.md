---
title: "Determinanty wskaźnika samobójstw: Przestrzenna analiza ekonometryczna wydatków na zdrowie, oczekiwanej długości życia i PKB per capita"
author: "Piotr Łukowski, Tomasz Kotliński, Jakub Gużewski"
date: "`r Sys.Date()`"
output:
  pdf_document:
  latex_engine: xelatex
toc: true
number_sections: true
---

```{r setup, include=FALSE}
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(car)
library(lmtest)
library(tseries)
library(viridis)
library(tidyr)
library(purrr)
library(readr)
library(whitestrap)
library(spatialreg)
```

# Lista państw i przypisanie kodów NUTS

```{r}
panstwa <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Czechia", "Denmark", "Estonia",
             "Finland", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy",
             "Latvia", "Lithuania", "Luxembourg", "Netherlands", "Norway", "Poland",
             "Portugal", "Slovak Republic", "Slovenia", "Spain", "Sweden", "Switzerland")

nazwa_do_nuts <- c(
  "Austria" = "AT", "Belgium" = "BE", "Bulgaria" = "BG", "Croatia" = "HR",
  "Czechia" = "CZ", "Denmark" = "DK", "Estonia" = "EE", "Finland" = "FI",
  "France" = "FR", "Germany" = "DE", "Greece" = "EL", "Hungary" = "HU",
  "Ireland" = "IE", "Italy" = "IT", "Latvia" = "LV", "Lithuania" = "LT",
  "Luxembourg" = "LU", "Netherlands" = "NL", "Norway" = "NO", "Poland" = "PL",
  "Portugal" = "PT", "Slovak Republic" = "SK", "Slovenia" = "SI",
  "Spain" = "ES", "Sweden" = "SE", "Switzerland" = "CH"
)
```

# Wczytanie mapy i funkcja pomocnicza

```{r}
eu_all <- st_read("C:/Users/gozia/Downloads/OneDrive_1_10.06.2025/NUTS_RG_01M_2021_4326.shp")

wczytaj_dane <- function(plik, nazwa_zmiennej, rok) {
  dane <- read.csv(plik, skip = 3, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(dane) <- gsub(" \\[YR[0-9]{4}\\]", "", colnames(dane))
  rok_str <- as.character(rok)
  dane <- dane[dane$`Country Name` %in% panstwa, ]
  if (!(rok_str %in% colnames(dane))) {
    stop(paste("Brak danych dla roku", rok_str, "w pliku", plik))
  }
  dane <- dane %>% select(`Country Name`, all_of(rok_str)) %>%
    rename(!!nazwa_zmiennej := all_of(rok_str))
  return(dane)
}
```

# Funkcja analizy

```{r}
analiza_rok <- function(rok) {
  d1 <- wczytaj_dane("C:/Users/gozia/Desktop/API_SH.XPD.CHEX.PC.CD_DS2_en_csv_v2_4535.csv", "HEX", rok)
  d2 <- wczytaj_dane("C:/Users/gozia/Desktop/API_SP.DYN.LE00.IN_DS2_en_csv_v2_2579.csv", "LEX", rok)
  d3 <- wczytaj_dane("C:/Users/gozia/Desktop/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2702.csv", "PKBC", rok)
  d4 <- wczytaj_dane("C:/Users/gozia/Desktop/API_SH.STA.SUIC.P5_DS2_en_csv_v2_7536.csv", "SCD", rok)
  
  dane_all <- reduce(list(d1, d2, d3, d4), full_join, by = "Country Name") %>%
    mutate(NUTS_ID = nazwa_do_nuts[`Country Name`])
  
  eu <- eu_all %>%
    filter(NUTS_ID %in% dane_all$NUTS_ID) %>%
    left_join(dane_all, by = "NUTS_ID") %>%
    mutate(
      log_HEX = log(HEX),
      log_LEX = log(LEX),
      log_SCD = log(SCD),
      log_PKBC = log(PKBC)
    ) %>%
    filter(complete.cases(log_HEX, log_LEX, log_PKBC, log_SCD))
  
  nb_q1 <- poly2nb(eu, queen = TRUE)
  islands <- which(card(nb_q1) == 0)
  if (length(islands) > 0) {
    eu <- eu[-islands, ]
    nb_q1 <- poly2nb(eu, queen = TRUE)
  }
  Wq1 <- nb2listw(nb_q1, style = "W", zero.policy = TRUE)
  
  model <- log_SCD ~ log_LEX + log_PKBC + log_HEX
  OLS <- lm(model, data = eu)
  SEM_model <- errorsarlm(model, data = eu, listw = Wq1)
  
  res <- residuals(OLS)
  bptest_res <- bptest(OLS, studentize = FALSE)
  wt <- white_test(OLS)
  jb <- jarque.bera.test(res)
  shapiro <- shapiro.test(res)
  moran <- lm.morantest(OLS, listw = Wq1)
  rstests <- lm.RStests(OLS, listw = Wq1, test = "all")
  
  local_moran <- localmoran(res, listw = Wq1)
  eu$local_I <- local_moran[, 1]
  eu$pval_I <- local_moran[, 5]
  
  p <- ggplot(eu) +
    geom_sf(aes(fill = local_I)) +
    scale_fill_viridis_c(option = "plasma", name = "Local Moran's I") +
    labs(title = paste("Lokalny I Morana (reszty OLS), rok", rok),
         caption = "Obliczono na podstawie reszt modelu OLS") +
    theme_minimal()
  print(p)
  
  cat("\n\n========== ANALIZA DLA ROKU:", rok, "==========\n")
  cat("\n--- OLS ---\n"); print(summary(OLS))
  cat("\n--- Testy diagnostyczne OLS ---\n")
  cat("Breusch-Pagan:\n"); print(bptest_res)
  cat("White:\n"); print(wt)
  cat("Jarque-Bera:\n"); print(jb)
  cat("Shapiro-Wilk:\n"); print(shapiro)
  cat("\n--- Testy przestrzenne ---\n")
  cat("Moran:\n"); print(moran)
  cat("LM tests:\n"); print(rstests)
  cat("\n--- SEM ---\n"); print(summary(SEM_model))
  cat("\n--- Log-Likelihood ---\n")
  cat("OLS:", as.numeric(logLik(OLS)), " | SEM:", as.numeric(logLik(SEM_model)), "\n")
  cat("\n--- AIC ---\n")
  print(AIC(OLS, SEM_model))
  
  return(list(
    rok = rok,
    OLS = OLS,
    diagnostyka = list(
      bptest = bptest_res,
      white = wt,
      jarque_bera = jb,
      shapiro = shapiro
    ),
    testy_przestrzenne = list(
      moran = moran,
      rstests = rstests
    ),
    SEM = SEM_model,
    logLik = c(OLS = as.numeric(logLik(OLS)), SEM = as.numeric(logLik(SEM_model))),
    AIC = AIC(OLS, SEM_model)
  ))
}
```

# Uruchomienie analizy dla kilku lat

```{r}
lata <- c(2016, 2018, 2020)
wyniki <- lapply(lata, analiza_rok)
```

# Tabela porównawcza

```{r}
tabela_porownawcza <- do.call(rbind, lapply(wyniki, function(x) {
  data.frame(
    Rok = x$rok,
    AIC_OLS = x$AIC[1],
    AIC_SEM = x$AIC[2],
    LogLik_OLS = x$logLik["OLS"],
    LogLik_SEM = x$logLik["SEM"]
  )
}))
print(tabela_porownawcza)
```
# Wstęp – Charakterystyka danych i zakres analizy
Celem niniejszego badania była identyfikacja oraz ilościowa ocena determinant wpływających na poziom samobójstw w krajach europejskich w ujęciu przestrzennym. Analiza miała na celu nie tylko wykrycie zależności statystycznych pomiędzy wybranymi zmiennymi społeczno-ekonomicznymi a wskaźnikiem samobójstw, ale również uchwycenie geograficznego wzorca ich występowania.

Uwzględnienie tła przestrzennego pozwala zaobserwować, czy i w jakich regionach Europy problem ten jest szczególnie nasilony, co może świadczyć o różnym poziomie świadomości społecznej, systemowego podejścia do zdrowia psychicznego lub odmiennym dostępie do profilaktyki i opieki psychiatrycznej.

Dzięki zastosowaniu narzędzi ekonometrii przestrzennej możliwe było określenie, czy kraje sąsiadujące cechują się zbliżonym poziomem wskaźnika samobójstw, oraz czy wpływ determinant (np. długości życia, wydatków na zdrowie, zamożności) układa się w określone wzory geograficzne. Wyniki mogą wskazywać nie tylko na ekonomiczne i zdrowotne uwarunkowania zjawiska, ale również pośrednio na świadomość problemu i priorytety polityki publicznej w poszczególnych krajach.

W szczególności badanie koncentruje się na wpływie trzech zmiennych: oczekiwanej długości życia (Life Expectancy – LEX), wydatków na zdrowie per capita (Health Expenditure – HEX) oraz produktu krajowego brutto per capita (GDP per capita – PKBC). Analizowana zmienna zależna to liczba samobójstw na 100 000 mieszkańców (Suicide Death Rate – SCD), przy czym wszystkie zmienne zostały poddane transformacji logarytmicznej.

Długość życia została potraktowana jako wskaźnik ogólnego stanu zdrowia społeczeństwa, odzwierciedlający jakość opieki medycznej, styl życia i skuteczność polityk zdrowotnych. Wydatki na zdrowie per capita interpretowano jako miarę dostępności i jakości systemu opieki zdrowotnej, który może mieć istotne znaczenie w kontekście profilaktyki i leczenia zaburzeń psychicznych. Natomiast PKB per capita uznano za wskaźnik poziomu zamożności społeczeństwa, który potencjalnie wpływa na dobrostan psychiczny poprzez warunki życia, bezpieczeństwo socjalne i dostęp do usług wspierających zdrowie psychiczne.
Do badania wykorzystano dane panelowe pochodzące z otwartych zasobów Banku Światowego oraz zbioru danych geograficznych w formacie NUTS (Nomenclature of Territorial Units for Statistics), które umożliwiły analizę przestrzenną z zastosowaniem modeli ekonometrycznych. Dane zostały przefiltrowane do poziomu krajowego (NUTS 0) i obejmują 26 państw europejskich, m.in. Niemcy, Francję, Włochy, Polskę, Hiszpanię, Szwecję, Norwegię, Czechy, Szwajcarię, Austrię, Grecję i inne.

Analizę przeprowadzono dla trzech punktów czasowych: 2016, 2018 oraz 2020 roku. Wybór tych lat podyktowany był dostępnością kompletnych danych dla wszystkich analizowanych zmiennych, a także chęcią uchwycenia ewentualnych zmian w czasie. Każdorazowo dane z danego roku zostały zintegrowane z warstwą przestrzenną, a następnie wykorzystane do estymacji modeli klasycznych (OLS) oraz przestrzennych (SEM – Spatial Error Model), co umożliwiło dokładniejsze ujęcie przestrzennej zależności między krajami.

# Uzasadnienie wyboru zjawiska oraz determinant
Zdecydowaliśmy się na analizę przestrzennego zróżnicowania wskaźnika samobójstw (Suicide Death Rate – SCD) w krajach europejskich, ponieważ jest to temat niezwykle istotny społecznie, szczególnie w dzisiejszych czasach. W ostatnich latach obserwuje się rosnącą liczbę debat publicznych i medialnych na temat zdrowia psychicznego oraz przyczyn samobójstw. Zjawisko to coraz częściej staje się przedmiotem refleksji społecznej, a jego zrozumienie – zarówno przez specjalistów, jak i obywateli – jest kluczowe dla tworzenia skutecznych strategii profilaktycznych i wsparcia.

W kontekście tej rosnącej świadomości społecznej uznaliśmy, że temat ten zasługuje na szczególną uwagę również w ujęciu naukowym i statystycznym. Postawiliśmy sobie pytanie: jakie czynniki ekonomiczne i zdrowotne mogą wpływać na decyzję jednostki o odebraniu sobie życia? Choć przyczyny samobójstw są zawsze złożone i wielowymiarowe, uznaliśmy, że trzy konkretne zmienne – oczekiwana długość życia (LEX), wydatki na zdrowie per capita (HEX) oraz PKB per capita (PKBC) mogą odgrywać istotną rolę w kształtowaniu poziomu samobójstw w poszczególnych krajach.

# Opis metodyki badawczej
W celu oszacowania zależności między wskaźnikiem samobójstw a wybranymi zmiennymi społeczno-ekonomicznymi, zastosowano podejście modelowania ekonometrycznego, uzupełnione o analizę przestrzenną. Podstawowym modelem była klasyczna regresja liniowa (OLS), którą następnie porównano z modelem przestrzennym typu SEM (Spatial Error Model), uwzględniającym autokorelację w składniku losowym.

$$
\log(\text{SCD}_i) = \beta_0 + \beta_1 \log(\text{LEX}_i) + \beta_2 \log(\text{PKBC}_i) + \beta_3 \log(\text{HEX}_i) + \varepsilon_i
$$


Po estymacji modelu OLS zastosowano szereg testów weryfikujących założenia klasycznej regresji liniowej:

Test Breuscha-Pagana oraz test White’a – oceniają obecność heteroskedastyczności (czyli niejednorodności wariancji reszt),

Test Jarque-Bera oraz test Shapiro-Wilka – badają normalność rozkładu reszt,

Test Morana (Global Moran’s I) – sprawdza przestrzenną autokorelację reszt.

Występowanie istotnej przestrzennej autokorelacji stanowi przesłankę do zastosowania modeli przestrzennych, ponieważ narusza ono podstawowe założenie niezależności błędów w modelu OLS.

Model przestrzenny SEM (Spatial Error Model):
Ze względu na wykrycie istotnej przestrzennej autokorelacji reszt we wszystkich analizowanych latach, zastosowano model przestrzenny typu SEM, który uwzględnia strukturę zależności przestrzennych w błędach. Model SEM ma postać:

$$
\log(\text{SCD}_i) = \beta_0 + \beta_1 \log(\text{LEX}_i) + \beta_2 \log(\text{PKBC}_i) + \beta_3 \log(\text{HEX}_i) + u_i
$$

$$
u_i = \lambda W u + \varepsilon_i
$$

#  Wyniki estymacji modeli OLS i SEM
W ramach analizy przeprowadzono estymację dwóch typów modeli: klasycznego modelu regresji liniowej (OLS) oraz modelu przestrzennego typu SEM (Spatial Error Model). Celem było sprawdzenie, jak trzy wybrane zmienne – oczekiwana długość życia, wydatki na zdrowie oraz PKB per capita – wpływają na poziom samobójstw w krajach europejskich. Analiza została przeprowadzona dla trzech lat: 2016, 2018 i 2020.

We wszystkich analizowanych latach najbardziej stabilną i istotną statystycznie zmienną okazała się oczekiwana długość życia. Wyniki wskazują, że w krajach, gdzie ludzie żyją dłużej, odnotowuje się niższe wskaźniki samobójstw. Taka zależność była widoczna zarówno w modelach klasycznych, jak i przestrzennych. Wskazuje to na ogólną zależność między lepszym stanem zdrowia publicznego a niższym poziomem zachowań samobójczych.

Zmienna dotycząca wydatków na zdrowie nie była istotna w modelach OLS, ale zyskiwała znaczenie w modelach przestrzennych SEM. Oznacza to, że po uwzględnieniu zależności geograficznych (czyli faktu, że kraje sąsiednie mogą na siebie wpływać), większe nakłady na zdrowie mogą mieć pozytywny wpływ na ograniczenie samobójstw – choć relacja ta nie jest jednoznaczna i może wynikać również z lepszej diagnostyki problemów psychicznych w bardziej rozwiniętych systemach opieki zdrowotnej.

Z kolei PKB per capita, czyli poziom zamożności danego kraju, nie wykazał istotnego wpływu na poziom samobójstw w żadnym z modeli. Może to oznaczać, że sama sytuacja materialna nie jest wystarczającym wyjaśnieniem zróżnicowania wskaźników samobójstw między państwami europejskimi.

Co ważne, modele przestrzenne SEM okazały się lepiej dopasowane do danych niż modele OLS. W każdym z analizowanych lat modele SEM wykazały lepsze wskaźniki jakości dopasowania, co sugeruje, że uwzględnienie zależności przestrzennych – czyli faktu, że kraje sąsiednie mogą mieć podobne cechy – jest uzasadnione i pozwala lepiej zrozumieć badane zjawisko.

# Diagnostyka modelu i analiza przestrzenna
Po oszacowaniu modeli OLS i SEM przeprowadzono szereg testów diagnostycznych, które miały na celu ocenę poprawności przyjętych założeń oraz sprawdzenie, czy uwzględnienie zależności przestrzennych było uzasadnione. Skupiono się na analizie reszt modelu – czyli różnic pomiędzy wartościami rzeczywistymi a przewidywanymi przez model.

Testy klasyczne (dla modelu OLS)
W pierwszej kolejności zastosowano testy oceniające jakość modelu OLS. Testy Breuscha-Pagana oraz White’a pozwoliły sprawdzić, czy reszty modelu mają jednakową zmienność (czyli czy nie występuje tzw. heteroskedastyczność). W większości przypadków nie stwierdzono problemów w tym zakresie – wyjątkiem był rok 2020, gdzie test White’a wykazał, że wariancja reszt nie jest stała, co może wpływać na wiarygodność wyników.

Kolejnym krokiem była ocena normalności rozkładu reszt – czyli sprawdzenie, czy rozkład błędów jest zbliżony do rozkładu normalnego, co jest jednym z założeń klasycznej regresji. Zarówno test Shapiro-Wilka, jak i Jarque-Bera wskazały, że reszty są rozłożone prawidłowo, a więc modele OLS nie naruszają tego założenia.

Testy przestrzenne – wykrycie zależności geograficznych
Najważniejszą częścią diagnostyki było jednak sprawdzenie, czy występuje autokorelacja przestrzenna – czyli sytuacja, w której wyniki w jednym kraju są podobne do wyników w krajach sąsiednich. Do tego celu zastosowano test Morana (Global Moran’s I), który we wszystkich trzech analizowanych latach wykazał istotne statystycznie zależności przestrzenne. Oznacza to, że klasyczny model OLS pomijał ważne wzorce przestrzenne w danych.

Dodatkowo zastosowano tzw. testy LM (Lagrange Multiplier), które pozwalają określić, jaki typ modelu przestrzennego będzie najodpowiedniejszy. Wyniki tych testów wskazywały na przewagę modelu SEM – czyli modelu z przestrzenną autokorelacją błędu – nad modelem z efektem lagowym (gdzie wpływ przestrzenny występuje bezpośrednio w wartości zmiennej zależnej).

Potwierdza to, że błędy modelu klasycznego są w dużej mierze powiązane z lokalizacją geograficzną – np. kraje położone obok siebie wykazują zbliżone odchylenia od wartości oczekiwanych. Pominięcie tego zjawiska prowadziłoby do zaniżenia wiarygodności wyników modelu O


# Wnioski i implikacje

Przeprowadzona analiza miała na celu określenie, które czynniki społeczno-ekonomiczne determinują poziom samobójstw w krajach europejskich oraz czy występują przestrzenne zależności w tym zjawisku. W badaniu uwzględniono trzy kluczowe zmienne niezależne: oczekiwaną długość życia (LEX), wydatki na zdrowie per capita (HEX) oraz PKB per capita (PKBC). W analizie wykorzystano dane dla 26 państw europejskich z trzech wybranych lat: 2016, 2018 i 2020. Zastosowano klasyczny model regresji liniowej (OLS) oraz przestrzenny model błędu (SEM), umożliwiający uwzględnienie autokorelacji przestrzennej w resztach modelu.

Wyniki estymacji jednoznacznie wskazują, że najsilniejszym i najbardziej stabilnym predyktorem wskaźnika samobójstw jest oczekiwana długość życia. We wszystkich analizowanych latach zmienna ta była istotna statystycznie, a jej wpływ negatywny, co oznacza, że w krajach, gdzie mieszkańcy żyją dłużej, odnotowuje się niższe poziomy samobójstw. Może to wynikać z wyższego ogólnego poziomu zdrowia publicznego, lepszej jakości życia, większego dostępu do systemu opieki zdrowotnej i wsparcia psychicznego.

Zmienna dotycząca wydatków na zdrowie per capita nie była istotna w modelach OLS, natomiast uzyskała istotność w modelach przestrzennych SEM. Świadczy to o tym, że relacja między nakładami na ochronę zdrowia a samobójstwami może być uzależniona od lokalnego kontekstu geograficznego. W krajach sąsiadujących, gdzie występuje podobna infrastruktura zdrowotna lub podobne wyzwania demograficzne, efekty wydatków na zdrowie mogą się przestrzennie nakładać. Może to także oznaczać, że większe wydatki wiążą się z lepszą diagnostyką i zgłaszalnością przypadków, co wpływa na sposób raportowania danych.

W przypadku PKB per capita nie stwierdzono istotności w żadnym z modeli. Wskazuje to, że poziom zamożności nie tłumaczy bezpośrednio różnic w poziomie samobójstw między krajami europejskimi. Sugeruje to potrzebę większego uwzględnienia czynników psychologicznych, społecznych oraz kulturowych w badaniach nad tym zjawiskiem.

Kluczowym elementem analizy było zastosowanie narzędzi ekonometrii przestrzennej. Testy diagnostyczne wykazały istotną autokorelację przestrzenną reszt w modelach OLS, co potwierdziły test Morana oraz testy Lagrange’a. To oznacza, że klasyczny model regresji pomijał istotne wzorce przestrzenne w danych, przez co jego wyniki mogły być zniekształcone. Modele przestrzenne SEM, które uwzględniają zależności przestrzenne w składniku losowym, okazały się trafniejsze – wykazywały niższe wartości AIC oraz wyższy poziom logarytmu wiarygodności, co potwierdza ich lepsze dopasowanie do danych.

Dodatkowym, ważnym elementem analizy była wizualizacja lokalnego wskaźnika Morana (Local Moran’s I) dla reszt modelu OLS. W każdej z trzech analizowanych lat (2016, 2018, 2020) zaobserwowano wyraźne lokalne skupisko o wysokim poziomie autokorelacji przestrzennej w rejonie Europy Południowo-Wschodniej, szczególnie w okolicach Grecji i Bałkanów. Region ten systematycznie odchylał się od wartości przewidywanych przez model – był obszarem o wysokich, istotnych lokalnych wartościach I Morana. Może to świadczyć o istnieniu lokalnych czynników kulturowych, społecznych lub instytucjonalnych, które nie zostały ujęte w modelu, a które istotnie wpływają na poziom samobójstw. Powtarzalność tego wzorca w czasie sugeruje, że nie mamy do czynienia z przypadkową anomalią, lecz z trwałą zależnością przestrzenną.
