package zapoctovy_program;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

// TODO udělat maven project
// TODO udělat testy na co nejvíc metod


public class Program {
    public static void main(String[] args){
        int gapPentalty = 0; String function = "sw";
        Path refSeqPath = Path.of(""); Path compSeqPath = Path.of("-");

        // TODO Parsování command line


        // Parsování matice
        // TODO rozšířit default matice o extra písmena
        // TODO opravit user-defined matrixIndexer
        MatrixReader matrixReader = new MatrixReader("","");

        // Parsování FASTA souborů
        // TODO dokumentace
        FastaParse fP = new FastaParse(refSeqPath,compSeqPath);

        // Alignment
        // TODO dokumentace
        List<Alignment> alignments = new ArrayList<>();
        for (Sequence seq:fP.compareSeq) {
            alignments.add(new Alignment(fP.refSeq,seq,gapPentalty,matrixReader,function));
        }

        // Output
        Output.toFile(alignments,80);
    }
}


class FastaParse {
    Sequence refSeq;
    Sequence[] compareSeq;

    /**
     *
     * @param refSeqPath
     * @param compSeqPath
     */
    public FastaParse(Path refSeqPath, Path compSeqPath){
        refSeq = parseFASTA(refSeqPath)[0]; // kdyby jich tam bylo více, tak se vezme jen tu první
        compareSeq = parseFASTA(compSeqPath);
    }

    /**
     *
     * @param path
     * @return
     */
    Sequence[] parseFASTA(Path path) {
        List<Sequence> seqList = new ArrayList<Sequence>();
        // Procházíme soubor a hledáme jednotlivé sekvence
        try (Scanner sc = new Scanner(path);){
            Sequence currSeq = null;
            while(sc.hasNext()){
                String line = sc.nextLine();
                // začátek nové sekvence
                if(line.charAt(0) == '>'){
                    if(currSeq != null){
                        // když narazíme na novou sekvenci, tak tu starou uzavřeme a uložíme
                        currSeq.end();
                        seqList.add(currSeq);
                    }
                    // definovní nové sekvence a headeru
                    currSeq = new Sequence(line);
                }
                // přidávání řádků sekvence
                if(currSeq != null){
                    currSeq.seqList.add(line);
                }
            }
        } catch (Exception e) {
            System.out.println("File not found at: " + path);
            System.exit(1);
        }
        if(seqList.isEmpty()){
            // konec, pokud se nanajde žádná sekvence v souboru
            System.out.println("No sequence found in: " + path);
            System.exit(1);
        }

        return seqList.toArray(new Sequence[0]);
    }
}

class Output{
    /**
     * Vypíše výsledky všech alignmentů do specifikovaného souboru. <br>
     * Pokud již existuje soubor se stejným jménem, tak má uživatel možnost starý přepsat nebo
     * si vybrat jiné jméno souboru <br>
     * Vždy nejdříve vypíše oba headery sekvencí a pak postupně vypisuje alignment.
     *
     * @param alignments list všech alignmentů
     * @param path cesta k souboru, kam se vypíší výsledky
     * @param lineLen délka jedné řádky alignmentu v souboru
     */
    static void toFile(List<Alignment> alignments, Path path, int lineLen){
        // když je specifikovaný path uživatelem
        File file = path.toFile();

        // vytvoření souboru
        try { boolean created = file.createNewFile();
            /*while(!created) {
                // loop, který zajistí, že se vytvoří nový output soubor
                // TODO asi nebude fungovat, protože přidáváme na konec jména, tedy na extention, ne do jména souboru
                file = new File(file.toPath() + "-");
                created = file.createNewFile();
            }
            System.out.println("Output file already exists, new file was created as " + file);*/
            if(!created){
                System.out.println("Output file already exists, type y/n whether you want to delete the existing" +
                        "one or abort to choose a different name");
                String choice = System.console().readLine();
                if(choice.equals("y")){
                    file.delete();
                    file.createNewFile();
                }
                else{
                    System.exit(0);
                }
            }

        } catch (Exception e){
            System.out.println("Output file error");
            System.exit(1);
        }

        // vložení do souboru alignemnt
        try {
            FileWriter fileWriter = new FileWriter(file);
            for(Alignment a: alignments){
                // zapsání headerů pro jasné poznání o které alignmenty se jedná
                fileWriter.write(a.seq1.header);
                fileWriter.write(a.seq2.header);
                fileWriter.write("\n");
                // zapsání alignmentu
                for(Tuple<String,String> t: a.toStringArrays(lineLen)){
                    fileWriter.write(t.item1);
                    fileWriter.write(t.item2);
                    fileWriter.write("\n");
                }
                // CIGAR a sckóre alignmentu
                fileWriter.write("CIGAR: " + a.cigarToStringArray(lineLen));
                fileWriter.write("Score: " + a.score);
                fileWriter.write("\n");
            }
            fileWriter.close();
        } catch (IOException e) {
            System.out.println("Error in writing to file: " + file );
            throw new RuntimeException(e);
        }

    }

    /**
     * Vypíše výsledky všech alignmentů do předdefinovaného souboru OUT.txt. <br>
     * Dále používá metodu toFile definovanou i s cestou k souboru. <br>
     * Je možné, že se stane, že již OUT.txt soubor exsituje, pak má uživatel možnost soubor
     * přepsat nebo si vybrat jiné, ne předdefinované, jméno.
     *
     * @param alignments list všech alignmentů
     * @param lineLen délka jedné řádky alignmentu v souboru
     */
    static void toFile(List<Alignment> alignments, int lineLen){
        Path defaultPath = Path.of(new File("").getAbsolutePath());
        defaultPath = Path.of(defaultPath + "OUT.txt");
        toFile(alignments,defaultPath,lineLen);
    }

    /**
     * Vypíše jenom headery, CIGAR string a skóre všech alignmentů do konzole.
     *
     * @param alignments list všech alignmentů
     */
    static void toCommandLine(List<Alignment> alignments){
        for(Alignment a:alignments){
            System.out.println(a.seq1.header);
            System.out.println(a.seq2.header);
            System.out.println("CIGAR: " + a.cigarToStringArray(80));
            System.out.println("Score: " + a.score);
            System.out.println();
        }
    }
}