package cz.cuni.mff.polakar;

import org.apache.commons.cli.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

// TODO udělat testy na co nejvíc metod
// TODO testFiles
// TODO     testfile s více sekvencemi
// TODO     ukázat, že to funguje na nějakých jednoduchých příkladech


public class Program {
    public static void main(String[] args){
        // zadefinování všech možných argumentů, které uživatel napíše do konzole
        int gapPenalty = -1; String alignmentFunction = null; String alignmentType = null;
        Path refSeqPath = null; Path compSeqPath = null;
        Path matrixPath = null; String matrixType = null;
        Path outPath = null; boolean toCommandLine = false;
        String MM = null;

        /*
        String defP = Paths.get("").toAbsolutePath() + "\\alignment\\";
        alignmentFunction = "nw"; alignmentType = "n";
        refSeqPath = Path.of(defP + "testFiles/nu1.fna"); compSeqPath = Path.of(defP + "testFiles/nu2.fna");
        //matrixType = "C:\\Users\\marti\\OneDrive\\Plocha\\ŠKOLA\\ŠKOLA 2023-2024\\java\\alignment\\src\\main\\java\\cz\\cuni\\mff\\polakar\\blosum62.txt";
        //matrixType = "\\alignment\\src\\main\\java\\cz\\cuni\\mff\\polakar\\blosum62.txt";
        //matrixType = Paths.get("").toAbsolutePath() + matrixType;
        matrixType = "defaultNucleotide2";
        matrixPath = Path.of(defP + "testFiles/myNuMatrix.txt");
        toCommandLine = false; */


        // command line options
        Options options = createOptions();
        CommandLineParser parser = new DefaultParser();

        try {
            CommandLine cmd = parser.parse(options, args);
            refSeqPath = Path.of(cmd.getOptionValue("refSeqFile"));
            compSeqPath = Path.of(cmd.getOptionValue("compSeqFile"));
            alignmentType = cmd.getOptionValue("alignmentType");
            alignmentFunction = cmd.getOptionValue("alignmentFunction");
            if(cmd.hasOption("matrixPath")){
                matrixPath = Path.of(cmd.getOptionValue("matrixPath"));
            }
            if(cmd.hasOption("matrixType")){
                matrixType = cmd.getOptionValue("matrixType");
            }
            if(cmd.hasOption("matrixMM")){
                MM = cmd.getOptionValue("matrixMM");
            }
            if(cmd.hasOption("outputFile")) {
                outPath = Path.of(cmd.getOptionValue("outputFile"));
            }
            if(cmd.hasOption("gapPenalty")){
                gapPenalty = Integer.parseInt(cmd.getOptionValue("gapPenalty"));
            }
            if(cmd.hasOption("cmd")){
                toCommandLine = true;
            }
        }
        catch (ParseException e) {
            System.err.println("Parsing failed. " + e.getMessage());
            System.exit(1);
        }

        // Parsování matice DONE
        // TODO otestovat MatrixReader(int match, int mismatch) přes commandline
        MatrixReader matrixReader = null;
        if(matrixPath == null){
            if(MM != null){
                String[] split = MM.split(",");
                matrixReader = new MatrixReader(Integer.parseInt(split[0]),Integer.parseInt(split[1]),alignmentType);
            }
            else if(matrixType == null){
                System.out.println("No matrix type was found");
                System.exit(1);
            }
            else{
                matrixReader = new MatrixReader(matrixType,alignmentType);
            }
        }
        else{
            matrixReader = new MatrixReader(matrixPath,alignmentType);
        }


        // Parsování FASTA souborů DONE
        FastaParse fP = new FastaParse(refSeqPath,compSeqPath);

        // Alignment DONE
        // TODO smazat pomocné funkce
        // TODO validate SW, zatim to moc nefunguje :(
        List<Alignment> alignments = new ArrayList<>();
        for (Sequence seq:fP.compareSeq) {
            alignments.add(new Alignment(fP.refSeq,seq,gapPenalty,matrixReader,alignmentFunction));
        }

        // Output DONE
        // TODO smazat smazat debugging
        if(outPath == null && !toCommandLine){
            Output.toFile(alignments,80);
        } else if (toCommandLine) {
            Output.toCommandLine(alignments,80);
        } else{
            Output.toFile(alignments,outPath,80);
        }
    }

    /** Vytvoření Options pro parosvání argumentů z konzole
     *
     * @return Options
     */
    static Options createOptions(){
        Options options = new Options();

        Option refSeqFile = Option.builder("r").option("r").longOpt("refSeqFile").argName("file")
                .desc("File of the reference sequence").required().hasArg()
                .build();
        Option compSeqFile = Option.builder("c").option("c").longOpt("compSeqFile").argName("file")
                .desc("File of the sequences that will be compared to the reference sequence").hasArg()
                .required().build();
        Option alignmentType = Option.builder("t").option("t").longOpt("alignmentType").argName("type")
                .desc("Type of alignment - protein/prot/p/aa, nucleotide/nu/n").hasArg()
                .required().build();
        Option alignmentFunction = Option.builder("f").option("f").longOpt("alignmentFunction").argName("name")
                .desc("Alignment function that will be used - sw (Smith-Waterman), nw (Needleman-Wunsch").hasArg()
                .required().build();
        Option matrixPath = Option.builder("p").option("p").longOpt("martixPath").argName("file")
                .desc("Path to custom matrix").hasArg()
                .required(false).build();
        Option matrix = Option.builder().option("m").longOpt("matrixType").argName("type").hasArg()
                .desc("Choice from the default matrices - defaultNucleotide (for nucleotide), blosum45/50/62/80/90 (for protein)")
                .required(false).build();
        Option matrixMM = Option.builder().option("i").longOpt("matrixMM").argName("int,int").hasArg()
                .desc("Define Match and Mismatch values for a nucleotide matrix")
                .required(false).build();
        Option outputFile = Option.builder().option("o").longOpt("outputFile").argName("file").hasArg()
                .desc("Path to the outputfile, not required")
                .required(false).build();
        Option gapPenalty = Option.builder().option("g").longOpt("gapPenalty").argName("int").hasArg()
                .desc("Gap penalty for the alignments")
                .required(false).build();
        Option toCommandLine = new Option("l","terminal",false,
                "Whether to write the output to the commandline");

        options.addOption(refSeqFile);
        options.addOption(compSeqFile);
        options.addOption(alignmentType);
        options.addOption(alignmentFunction);
        options.addOption(matrixPath);
        options.addOption(matrix);
        options.addOption(outputFile);
        options.addOption(gapPenalty);
        options.addOption(toCommandLine);

        return options;
    }
}


class FastaParse {
    Sequence refSeq;
    Sequence[] compareSeq;

    /** Třída pro parsování FASTA souborů se sekvencemi. <br>
     * Bere si na vstup oba vložené soubory a vrací sekvence
     *
     * @param refSeqPath cesta k souboru s referenční sekvencí
     * @param compSeqPath cesta k souboru s ostatními sekvencemi, které budou porovnávány s referenční
     */
    public FastaParse(Path refSeqPath, Path compSeqPath){
        refSeq = parseFASTA(refSeqPath)[0]; // kdyby jich tam bylo více, tak se vezme jen tu první
        compareSeq = parseFASTA(compSeqPath);
    }

    /** Parsování samotného FASTA souboru se sekvencemi. <br>
     * Pokud se v souboru žádná sekvence nenajde, tak program ukončí.
     *
     * @param path cesta k FASTA souboru
     * @return pole sekvencí, které se nachází v souboru
     */
    Sequence[] parseFASTA(Path path) {
        List<Sequence> seqList = new ArrayList<>();
        // Procházíme soubor a hledáme jednotlivé sekvence
        try (Scanner sc = new Scanner(path)){
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
                    continue;
                }
                // přidávání řádků sekvence
                if(currSeq != null){
                    currSeq.seqList.add(line);
                }
            }
            if(currSeq != null) {
                currSeq.end();
                seqList.add(currSeq);
            }
        } catch (IOException e) {
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
            if(!created){
                System.out.println("Output file already exists, type y/n whether you want to delete the existing " +
                        "one or abort to choose a different name");
                Scanner sc = new Scanner(System.in);
                String choice = sc.nextLine();
                if(choice.equals("y")){
                    file.delete();
                    file.createNewFile();
                    System.out.println("Old file deleted and replaced with an empty file");
                }
                else{
                    System.exit(0);
                }
            }

        } catch (Exception e){
            System.out.println("Output file error" + e);
            System.exit(1);
        }

        // vložení do souboru alignemnt
        try {
            FileWriter fileWriter = new FileWriter(file);
            for(Alignment a: alignments){
                // zapsání headerů pro jasné poznání o které alignmenty se jedná
                fileWriter.write(a.seq1.header + "\n");
                fileWriter.write(a.seq2.header + "\n");
                fileWriter.write("\n");
                // zapsání alignmentu
                for(Tuple<String,String> t: a.toStringArrays(lineLen)){
                    fileWriter.write(t.item1 + "\n");
                    fileWriter.write(t.item2 + "\n");
                    fileWriter.write("\n");
                }
                // CIGAR a sckóre alignmentu
                fileWriter.write("CIGAR string:" + "\n");
                for(String s:a.cigarToStringArray(lineLen)){
                    fileWriter.write(s + "\n");
                }
                fileWriter.write("Score: " + a.score + "\n");
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
        defaultPath = Path.of(defaultPath + "\\alignment\\OUT.txt");
        //System.out.println(defaultPath); TODO
        toFile(alignments,defaultPath,lineLen);

    }

    /**
     * Vypíše jenom headery, CIGAR string a skóre všech alignmentů do konzole.
     *
     * @param alignments list všech alignmentů
     */
    static void toCommandLine(List<Alignment> alignments, int lineLen){
        for(Alignment a:alignments){
            System.out.println(a.seq1.header);
            System.out.println(a.seq2.header);

            // TODO smazat, do konzole se nevypisují alignmenty
            for(Tuple<String,String> t: a.toStringArrays(80)){
                System.out.println(t.item1);
                System.out.println(t.item2);
                System.out.println();
            }
            System.out.println("CIGAR:");
            for(String s:a.cigarToStringArray(lineLen)){
                System.out.println(s);
            }
            System.out.println("Score: " + a.score);
            System.out.println();
        }
    }
}