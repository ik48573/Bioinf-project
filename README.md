# Bioinf-project

Preduvjeti
  - gcc 9

Nakon preuzimanja, kompajlira se kroz linux okruženje preko komandne linije naredbom:
  g++-9 FileReader.cpp -std=c++17 -Ispoa/include -Lbuild/lib/ -lspoa -o example
      

Pokretanje:

"./example J29_B_CE_IonXpress_005.fastq"
  -- ovaj primjer generira 4 konsenzusa
  
"./example J30_B_CE_IonXpress_006.fastq"
  -- ovaj primjer kod poziva generate_consensus() [linija 334] javlja se greška "Segmentation fault(core dumped)"
  
