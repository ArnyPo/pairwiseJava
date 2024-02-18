package cz.cuni.mff.polakar;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class MatrixReader {
    Path matrixPath;
    String matrixType;
    char alignmentType;
    HashMap<Character,Integer> matrixIndexer = new HashMap<>();
    int[][] matrix;

    /**
     * Konstruktor pro třídu MatrixReader pro předdefinované matice
     *
     * @param matrixType typ matice (BLOSUM62 atp.)
     * @param alignmentType typ alignmentu - protein/nucleotide
     */
    public MatrixReader(String matrixType, String alignmentType){
        // jen pro předdefinované matice
        matrixPath = Paths.get(Paths.get("").toAbsolutePath().toString(), "\\alignment\\defaultMatrices\\", matrixType + ".txt");
        this.matrixType = matrixType;
        setAlignmentType(alignmentType);
        setMatrixIndexerDefault();
        setMatrix();

    }

    /**
     * Konstruktor pro třídu MatrixReader pro matici definovanou uživatelem
     *
     * @param matrixPath cesta k matici
     * @param alignmentType typ alignmentu - protein/nucleotide
     */
    public MatrixReader(Path matrixPath, String alignmentType){
        // pro user-defined matice
        this.matrixPath = matrixPath;
        setAlignmentType(alignmentType);
        setMatrixIndexerCustom();
        setMatrix();
    }

    /**
     * Konstruktor pro třídu MatrixReader pro vytvoření matice pro nucleotide alignment
     * zadává se jen match a mismatch a podle toho se vytvoří matice
     * @param match hodnota pro match mukleotidů
     * @param mismatch hodnota pro mismatch nukleotidů
     */
    public MatrixReader(int match, int mismatch, String alignmentType){
        setAlignmentType(alignmentType);
        setMatrixIndexerDefault();
        if(this.alignmentType != 'n'){
            System.out.println("Input alignment type is not compatible with using custom match/mismatch values");
            System.exit(1);
        }

        // vyplnime matici, na digonále jsou match, na ostatních jsou mismatch
        for(int i=0;i<4;i++){
            matrix[i] = new int[matrixIndexer.size()];
            for(int j=0;j<4;j++){
                if((i-j)==0){ // diagonála
                    matrix[i][j] = match;
                }
                else{
                    matrix[i][j] = mismatch;
                }
            }
        }
    }

    /**
     * Metoda ze souboru vytvoří matici formátu int[][] <br>
     * Když nenajde soubor, program končí
     */
    private void setMatrix(){
        try (Scanner scanner = new Scanner(matrixPath)){
            int i = 0;
            if(scanner.hasNext()) { // přeskočíme první řádek, ale kontrolujeme
                String line = scanner.nextLine();
                String[] header = Arrays.stream(line.split(" ")).filter(x -> !x.isEmpty()).toArray(String[]::new);
                if(header.length != matrixIndexer.size()){
                    System.out.println("Matrix type and alignment type are not compatible");
                    System.exit(1);
                }
            }
            while(scanner.hasNext()){
                String line = scanner.nextLine();
                String[] lineSplit = Arrays.stream(line.split(" ")).filter(x -> !x.isEmpty()).toArray(String[]::new);
                int[] arr = new int[lineSplit.length-1];
                for (int j=0; j<lineSplit.length-1; j++){
                    // první je vždycky písmeno, přidáváme 1
                    arr[j] = Integer.parseInt(lineSplit[j+1]);
                }
                matrix[i++] = arr;
            }
        }
        catch (IOException e){
            System.out.println("Matrix file not found at: " + matrixPath);
            System.exit(1);
        }
    }

    /**
     * Metoda hledá podle indexů znaků skóre v matici
     *
     * @param a,b znaky, pro které se hledá skóre
     * @return skóre z matice
     */
    public int reader(char a, char b){
        int indexA = matrixIndexer.get(a);
        int indexB = matrixIndexer.get(b);
        return matrix[indexA][indexB];
    }

    /**
     * Zkontroluje správnost aligmentType <br>
     * Pokud se proměnná neshoduje s žádnou očekávanou, tak ukončuje program. <br>
     * Do this.alignmentType zapisuje příslušný znak
     *
     * @param alignmentType typ alignmentu který zadává uživatel
     */
    private void setAlignmentType(String alignmentType) {
        if (new ArrayList<>(List.of("p","prot","protein","aa")).contains(alignmentType)) {
            this.alignmentType = 'p';
        }
        else if(new ArrayList<>(List.of("n","nu","nucleotide")).contains(alignmentType)){
            this.alignmentType = 'n';
        }
        else {
            System.out.println("Wrong alignment type, choose from nucleotide or protein");
            System.exit(1);
        }
    }

    /**
     * Nastaví matrixIndexer pro default matic
     */
    private void setMatrixIndexerDefault(){
        // pro default matice
        char[] nHeader = "ACGT".toCharArray();
        char[] pHeader = "ARNDCQEGHILKMFPSTWYVBJZX*".toCharArray();
        if (alignmentType == 'n'){
            for(int i=0; i<nHeader.length; i++){
                matrixIndexer.put(nHeader[i],i);
            }
        }
        else{
            for(int i=0; i<pHeader.length; i++){
                matrixIndexer.put(pHeader[i],i);
            }
        }
        matrix = new int[matrixIndexer.size()][];
    }

    /**
     * Nastaví matrixIndexer pro matice definované uživatelem
     */
    private void setMatrixIndexerCustom(){
        // pro user-defined matice
        // minimum písmen, která musí být obsažena
        // a extra písmena, která nemusí být obsažena
        String minP = "ARNDCQEGHILKMFPSTWYVBJZX";
        String extraP = "OU*";
        String minN = "ACGT";
        String extraN = "UiRYKMSWBDHVN";

        // check správnosti headeru
        String header = "";
        try {
            Scanner scanner = new Scanner(matrixPath);
            header = scanner.nextLine();
        }
        catch (Exception e){
            System.out.println(e);
            System.exit(1);
        }

        String finalHeader;
        if(alignmentType == 'n'){
            finalHeader = checkMatrixHeader(minN, extraN, header);
        }
        else{
            finalHeader = checkMatrixHeader(minP, extraP, header);
        }

        // vytvoření indexeru
        // nejdříve se raw header rozdělí na znaky a pak se bez prázných zase zpátky sloučí do jednoho
        //String headerString = String.join("",Arrays.stream(header.split("")).filter(x -> !x.isEmpty()).toArray(String[]::new)); TODO smazat
        String headerString = Arrays.stream(header.split("")).filter(x -> !x.isEmpty()).collect(Collectors.joining(""));
        for(int i=0; i<finalHeader.length(); i++){
            // procházím default header (finalHeader) a pro každý znak najdu, kde se nachází v custom headeru
            int pos = headerString.indexOf(finalHeader.charAt(i));
            matrixIndexer.put(headerString.charAt(pos),pos);
            //System.exit(0);
        }

        // vytoření matice
        matrix = new int[matrixIndexer.size()][];
    }

    /**
     * Kontrola, že ve skórovací matici jsou všechny potřebná písmena
     * @param min minimum, co je potřeba
     * @param extra extra písmena
     * @param header header z matice
     * @return header, který má v sobě všechny náležité písmena
     */
    private String checkMatrixHeader(String min, String extra, String header) {
        String finalHeader = ""; boolean hasMinLetters = true; boolean hasExtraLetters = true;
        String[] toCheckMinN = min.split("");
        String[] toCheckExtraN = extra.split("");
        for(String s:toCheckMinN){
            if (!header.contains(s)) {
                hasMinLetters = false;
                break;
            }
        }
        for(String s:toCheckExtraN){
            if (!header.contains(s)) {
                hasExtraLetters = false;
                break;
            }
        }
        if(hasMinLetters && hasExtraLetters){
            finalHeader = min + extra;
        }
        else if(hasMinLetters){
            finalHeader = min;
        }
        else {
            System.out.println("User-defined matrix is not compatible");
            System.exit(1);
        }
        return finalHeader;
    }
}
