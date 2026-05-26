# clustering-gene-variants

Grupiranje varijanti sekvenciranja jelena korištenjem višestrukog poravnanja sekvenci (MSA) i pohlepnog klasteriranja (engl. greedy clustering).

## Ovisnosti

Python 3.8+ i EMBOSS alat `needle`:

```bash
pip install -r requirements.txt
sudo apt install emboss
```

Python paketi: `biopython`, `edlib`, `numpy`, `pyspoa`

## Struktura projekta

```
data/
  fastq/                    # ulazne FASTQ datoteke (jedan uzorak po datoteci)
  J29B_expected.fasta       # referentna baza alela J29B
  J30B_expected.fasta       # referentna baza alela J30B
results/                    # generira se automatski pri pokretanju
analysis_results.csv        # rezultati (generirano s analyze.py)
```

## Pokretanje

Projekt se pokreće u tri koraka:

```bash
# 1. Obraditi uzorke za svaku referentnu bazu
python main.py config_J29.txt
python main.py config_J30.txt

# 2. Skupljanje rezultate poravnanja u CSV
python analyze.py
```

`main.py` obrađuje jedan uzorak prema jednoj referentnoj bazi. Za obradu više uzoraka pokrenuti jednom po uzorku.

## Konfiguracija

Primjer `config.txt`:

```
reference_file=data/J29B_expected.fasta
fragments_file=data/fastq/J29_B_CE_IonXpress_005.fastq
gapopen=10
gapextend=0.5
fragment_length=296
```

| Parametar         | Opis                                                         |
|-------------------|--------------------------------------------------------------|
| `reference_file`  | FASTA datoteka s referentnim alelima                         |
| `fragments_file`  | FASTQ datoteka s očitanjima jednog uzorka                    |
| `gapopen`         | Kazna za otvaranje praznine u needle poravnanju              |
| `gapextend`       | Kazna za produljenje praznine u needle poravnanju            |
| `fragment_length` | Duljina očitanja koja se uzimaju u obzir (za jelena: 296 bp) |

Uzorci s manje od 100 očitanja ciljane duljine automatski se isključuju. Klasteri s manje od 10 očitanja odbacuju se.

## Izlaz

| Datoteka                                    | Sadržaj                                              |
|---------------------------------------------|------------------------------------------------------|
| `results/*_representatives.fasta`          | Reprezentativne sekvence svih klastera jednog uzorka |
| `results/*_cluster*_representative.fasta`  | Reprezentativna sekvenca pojedinog klastera          |
| `results/*_cluster*_alignment.txt`         | needle poravnanje s referentnim alelima              |
| `analysis_results.csv`                      | Svi klasteri s najboljim alelom i % identičnosti     |
