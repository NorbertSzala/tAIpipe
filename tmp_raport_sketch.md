chciałbym się bardziej skupić na tym i jak to sensownie zaprojektować, ponieważ chcemy wziąć elementy z tej analizy i zrobić z tego analize dla dowolnego organizmu, to jest będziemy mieli różne królestwa dokmeny (czy to się przewija później w analizach i dlaczego???)) 



i teraz jaki informacje mogą byc interesujące jeżęli chodzi o ten indeks tAI i dlaczego (ponieważ chce zrozumieć co jest ciekawego w badaniu takiego zjawiska) 


# Sketch 

ponieważ na pewno widzę to tak, że:


## Wstęp / wstępne analizy QC ???? 

### Ogólny opis grup
CEL:
- upewnienie się co jest analizowane
METODOLOGIA:
- Spis Kohorty (Ingestion Summary): Dynamiczna tabela (pakiet DT w R) pokazująca liczbę próbek z podziałem na Domeny, Królestwa i Style Życia.

jakis ogólny opis tych grup które używamy itd. takie zestawienie, żeby wiedzieć na jakich danych działaliśmy, żeby osoba która dawała input była zapewniona co było wzięte pod uwagę w raporcie

### Analiza per genom
CEL:
- opisanie jakiś ogólnych statystykl per próbkę
- chcemy ocenić jakość pomiaru/danych takiej próbki na szybko
  - zwrócić przy okazji uwagi jeśli jakieś statystyki wyglądają skrajnie dziwnie
METODOLOGIA:
- agregujemy jakies informacje p[er próba ]
- Genomowe Wykresy Kontrolne (Sample QC): * Wykres rozrzutu (Scatter Plot): Wielkość genomu vs Całkowita liczba wykrytych genów tRNA. Służy do wykrywania anomalii (np. jeśli próbka ma za mało tRNA jak na swój rozmiar genomu, podświetlamy ją na czerwono jako ostrzeżenie).
  - Rozkład zawartości GC: Porównanie globalnego GC genomu z GC na trzeciej pozycji kodonu (GC3s) dla każdego gatunku, aby ocenić ogólny dryf mutacyjny.
MOŻE:
- jakby zrobić jakiś skrypt który wyciąga informacje z jakiś stron dotyczące tych gatunków (mamy grupy taksonomiczne itd.) i jakoś to opisuje w określony sposób 


***
Główne Analizy 
***

Tutaj byśmy sprawdzali jakies określone rzeczy, najchętniej  bbym chciał zrobić ta:
- w jaki sposób ty grupujesz te mezo itd.  (bo rozumiem że tutaj w innych skalach odpowiadamy to jest):
    - poziom białek itd.
    - poziom organizmu
    - poziom ogólnie jakieś większej domeny kohory??? \
- wtedy jakei analizy było by sensownei przeprowadzić:
  - najpierw bym chciał streszczenie tych które były i na jakie problemy odpowiadają, potem dać propozycje nowych

## Poziom Per drzewo taksonomiczne (no czyli per  grupy)_ 

## Poziom organizmu

## Poziom białka geny 

## Poziom kodonów tRNA 