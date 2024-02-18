package zapoctovy_program;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class MatrixReader {
    Path matrixPath;
    String matrixType;
    char alignmentType; // n nebo p
    HashMap<Character,Integer> matrixIndexer;
    int[][] matrix;

    /**
     * Konstruktor pro třídu MatrixReader pro předdefinované matice
     *
     * @param matrixType typ matice (BLOSUM62 atp.)
     * @param alignmentType typ alignmentu - protein/nucleotide
     */
    public MatrixReader(String matrixType, String alignmentType){
        // jen pro předdefinované matice
        matrixPath = Path.of("defaultMatrices" + matrixType);
        this.matrixType = matrixType;
        setAlignmentType(alignmentType);
        setMatrixIndexerDefault(); // TODO
        setMatrix();

    }

    /**
     * Konstruktor pro třídu MatrixReader pro matici definovanou živatelem
     *
     * @param matrixPath cesta k matici
     * @param matrixType ??SMAZAT?? typ matice (BLOSUM62 atp.)
     * @param alignmentType typ alignmentu - protein/nucleotide
     */
    public MatrixReader(Path matrixPath, String matrixType, String alignmentType){
        // pro user-defined matice
        // potřebuju matrixType??
        this.matrixPath = matrixPath;
        setAlignmentType(alignmentType);
        setMatrixIndexerCustom(); // TODO
        setMatrix();
    }

    /**
     * Metoda ze souboru vytvoří matici formátu int[][] <br>
     * Když nenajde soubor, program končí
     */
    void setMatrix(){
        try (Scanner scanner = new Scanner(matrixPath)){
            int i = 0;
            if(scanner.hasNext()){scanner.next();} // přeskočíme první řádku
            while(scanner.hasNext()){
                String line = scanner.nextLine();
                String[] lineSplit = line.split(" ");
                int[] arr = new int[lineSplit.length];
                for (int j=1; j<lineSplit.length; j++){
                    // první je vždycky písmeno, začínáme od 1
                    arr[j] = Integer.parseInt(lineSplit[j]);
                }
                matrix[i++] = arr;
            }
        }
        catch (Exception e){
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
    int reader(char a, char b){
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
    void setAlignmentType(String alignmentType) {
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
    void setMatrixIndexerDefault(){
        // pro default matice
        // TODO rozšířené matice o nestandartní písmena
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
    void setMatrixIndexerCustom(){
        // pro user-defined matice
        // minimum písmen, která musí být obsažena
        // a extra písmena, která nemusí být obsažena
        String minP = "ARNDCQEGHILKMFPSTWYVBJZX";
        String extraP = "OU*";
        String minN = "ACGT";
        String extraN = "UiRYKMSWBDHVN";

        // check správnosti headeru
        String header;
        try {
            Scanner scanner = new Scanner(matrixPath);
             header = scanner.nextLine();
        }
        catch (Exception e){
            throw new RuntimeException(e);
        }

        String finalHeader = "";
        boolean hasMinLetters = true;
        boolean hasExtraLetters = true;
        // TODO zjistit, jestli to dělá, to, co chci, aby to dělalo
        // prvaděpodobně ne, asi nebude správně fungovat .contains()
        if(alignmentType == 'n'){
            String[] toCheckMinN = minN.split(".");
            String[] toCheckExtraN = extraN.split(".");
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
            if(hasMinLetters){
                finalHeader = minN;
            }
            else if(hasMinLetters && hasExtraLetters){
                finalHeader = minN + extraN;
            }
            else {
                System.out.println("User-defined matrix is not compatible");
                System.exit(1);
            }
        }
        else{
            String[] toCheckMinP = minP.split(".");
            String[] toCheckExtraP = extraP.split(".");
            for(String s:toCheckMinP){
                if(!header.contains(s)){
                    hasMinLetters = false;
                }
            }
            for(String s:toCheckExtraP){
                if(!header.contains(s)){
                    hasExtraLetters = false;
                }
            }
            if(hasMinLetters){
                finalHeader = minP;
            }
            else if(hasMinLetters && hasExtraLetters){
                finalHeader = minP + extraP;
            }
            else {
                System.out.println("User-defined matrix is not compatible");
                System.exit(1);
            }
        }

        // vytvoření indexeru
        char[] charHeader = finalHeader.toCharArray();
        for(int i=0; i<header.length(); i++){
            matrixIndexer.put(charHeader[i],i);
        }

        // vytoření matice
        matrix = new int[matrixIndexer.size()][];
    }
}
