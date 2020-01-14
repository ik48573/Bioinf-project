# Bioinf-project

Upute za instalaciju i korištenje

Preduvjeti
  - gcc 9
	
Prije pokretanja, potrebno je u direktorij u kojem se nalazi FindConsensus.cpp datoteka dodati direktorij „fastq“ u kojem će se nalaziti datoteke fastq formata koje program obrađuje.

Nakon preuzimanja, kompajlira se kroz Linux okruženje preko komandne linije naredbom:
g++-9 FindConsensus.cpp -std=c++17 -Ispoa/include -Lbuild/lib/ -lspoa -o example

gdje –Ispoa/include predstavlja direktorij sa zaglavljima koje spoa koristi, -Lbuild/lib/ putanju do spoa biblioteke, a example proizvoljan naziv izvršne datoteke.

Primjer pokretanja: 
"./example J29_B_CE_IonXpress_005.fastq"
Tijekom izvođenja programa u terminalu ispisuje se pojedini korak tijeka izvođenja.

Nakon završetka programa u osnovnom direktoriju pojavljuje se direktorij „fasta“ koji u sebi sadrži fasta datoteke koje su rezultat izvođena programa. 

  
