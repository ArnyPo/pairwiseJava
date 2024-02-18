package zapoctovy_program;
import java.util.*;
public class Alignment {
    Sequence seq1;
    Sequence seq2;
    int gapPenalty;
    MatrixReader matrixReader;
    List<Tuple<Character,Character>> alignment = new LinkedList<Tuple<Character,Character>>();
    List<Tuple<Integer,Character>> cigarS = new LinkedList<Tuple<Integer,Character>>();
    int score;

    /**
     * Konstruktor pro třídu Alignment
     *
     * @param seq1 sekvence k porovnání
     * @param seq2 sekvence k porovnání
     * @param gapPenalty gap penalty pro obě funkce
     * @param matrixReader matice, která bude používána pro porovnávní
     * @param function funkce použitá pro alignment - Needleman-Wunsch (nw) nebo Smith-Watermann (sw)
     */
    public Alignment(Sequence seq1, Sequence seq2, int gapPenalty, MatrixReader matrixReader, String function){
        this.seq1 = seq1;
        this.seq2 = seq2;
        this.gapPenalty = gapPenalty;
        this.matrixReader = matrixReader;

        if(function.equals("nw")){
            needlemanWunsh();
        } else if (function.equals("sw")) {
            smithWaterman();
        }
        else {
            System.out.println("Wrong alignment function, function: \"" + function + "\" does not exist.");
            System.exit(1);
        }
    }

    /**
     * Needleman-Wunsch algoritmus pro lokální pairwise alignment <br>
     * f(i,j) = max {f(i-1,j-1)+m =A, f(i-1,j)+d =B, f(i,j-1)+d =C} <br>
     * d - gap penalty - výchozí hodnota -10 <br>
     * m(a,b) - similarity function <br>
     */
    void needlemanWunsh(){
        //if(gapPenalty==0){ gapPenalty=-10;}

        Match[][] plot = new Match[seq1.useLen][seq2.useLen];

        // inicializace
        plot[0][0] = new Match(0);
        for(int i=0; i< seq1.useLen; i++){
            plot[i][0] = new Match(i*gapPenalty,-1,0);
        }
        for(int i=0; i<seq2.useLen; i++){
            plot[0][i] = new Match(i*gapPenalty, 0,-1);
        }

        // vyplnění tabulky
        for(int i=1; i<seq1.useLen; i++){
            for(int j=1; j< seq2.useLen; j++){
                int M = matrixReader.reader(seq1.sequence.charAt(i-1),seq2.sequence.charAt(j-1));
                int A = plot[i-1][j-1].score + M;
                int B = plot[i-1][j].score + gapPenalty;
                int C = plot[i][j-1].score + gapPenalty;
                // výběr nejlepšího směru s preferencí
                plot[i][j] = new Match(bestDir(A,B,C,'A'));
            }
        }
        Tuple<Integer,Integer> start = new Tuple<Integer,Integer>(seq1.size, seq2.size);
        score = plot[seq1.useLen][seq2.useLen].score;
        align(plot,start);
    }

    /**
     * Smith-Waterman algoritmus pro globální alignment
     * f(i,j) = max{f(i-1,j-1) + m(a,b) =A, max{f(i-k,j)+d} =B, max{f(i,j-k)+d} =C, 0}
     * d - gap penalty
     * m(a,b) - similarity function
     */
    void smithWaterman(){
        //if(gapPenalty==0){ gapPenalty=-1;}
        Match[][] plot = new Match[seq1.useLen][seq2.useLen];
        Direction maxPoint = new Direction();

        // inicializace
        plot[0][0] = new Match(0);
        for(int i=0; i< seq1.useLen; i++){
            plot[i][0] = new Match(0);
        }
        for(int i=0; i<seq2.useLen; i++){
            plot[0][i] = new Match(0);
        }
        // vyplnění tabulky
        for(int i=1; i<seq1.useLen; i++){
            for(int j=1; j<seq2.useLen; j++){
                // jednotlivé části funkce
                int M = matrixReader.reader(seq1.sequence.charAt(i-1),seq2.sequence.charAt(j-1));
                int A = plot[i-1][j-1].score + M;
                int B = plot[i-1][j].score + j*gapPenalty;
                int C = plot[i][j-1].score + j*gapPenalty;
                // výběr nejlepšího směru s preferencí
                Direction bD = bestDir(A,B,C,'A');
                // byl bD gap?
                boolean gap = bD.getScore() == B || bD.getScore() == C;

                if(bD.getScore() >= 0){
                    if(bD.getScore() > maxPoint.getScore()){
                        maxPoint.modify(i, j, bD.getScore());
                    }
                    if(gap){
                        int gapLen = plot[i+bD.i][j+bD.j].gapLen + 1;
                        plot[i][j] = new Match(bD,gapLen);
                    }
                    else{
                        plot[i][j] = new Match(bD);
                    }
                }
                else{
                    plot[i][j] = new Match(0);
                }
            }
        }
        score = maxPoint.getScore();
        align(plot,maxPoint.getDirTup());
    }

    /**
     * Zde se provádí samotný alignment <br>
     * Prochází předpočítanou tabulku a podle instancí Direction se posouvá
     * a zapisuje příslušné hodnoty do alignment.
     * Zároveň vytváří CIGAR string.
     *
     * @param plot tabulka předpočítaných hodnot pro každý prvek sekvencí
     * @param start indexy, odkud má alignment začít
     */
    void align(Match[][] plot, Tuple<Integer,Integer> start){
        Match curr = plot[start.item1][start.item2];

        while(!curr.dir.equals(0,0)){
            Direction currDir = curr.dir;
            int di = start.item1 + currDir.i;
            int dj = start.item2 + currDir.j;
            curr = plot[di][dj];

            // podle toho, kam ukazuje currDir, tak vyplňujeme tabulku
            if(currDir.i == -1 && currDir.j == -1){
                alignment.addFirst(new Tuple<>(seq1.sequence.charAt(start.item1-1),
                                seq2.sequence.charAt(start.item2-1)));
                generateCIGAR('M');
            } else if (currDir.i == -1 && currDir.j == 0) {
                alignment.addFirst(new Tuple<>(seq1.sequence.charAt(start.item1-1),'-'));
                generateCIGAR('D');
            }
            else{
                alignment.addFirst(new Tuple<>('-',seq2.sequence.charAt(start.item2-1)));
                generateCIGAR('I');
            }
            start.update(di, dj);
        }
    }

    /**
     * Vrací instaci třídy Direction jako nejlepší možný směr, kam se může alignment vydat. <br>
     * Uživetel ale má možnost si zadefinovat preferenci směru, pokud by jich bylo více stejně dobrých.
     *
     * @param A hodnota pro match
     * @param B hodnota pro mismatch
     * @param C hodnota pro mismatch
     * @param preference preference hodnoty, která se zapíše
     * @return instance třídy Direction
     */
    Direction bestDir(int A, int B, int C, char preference){
        // preference A/B/C
        if(preference == 'A'){
            if(A>=B && A>=C) {return new Direction(-1,-1,A);}
            if(B>=A && B>=C) {return new Direction(-1,0,B);}
            return new Direction(0,-1,C);
        } else if (preference == 'B') {
            if(B>=A && B>=C) {return new Direction(-1,0,B);}
            if(A>=B && A>=C) {return new Direction(-1,-1,A);}
            return new Direction(0,-1,C);
        }
        else{
            if(C>=A && C>=B) {return new Direction(0,-1,C);}
            if(A>=B && A>=C) {return new Direction(-1,-1,A);}
            return new Direction(-1,0,B);
        }
    }

    /**
     * Metoda postupně generující CIGAR řetězec podle alignmentu
     *
     * @param type typ, který se má přidat do řetězce
     */
    void generateCIGAR(char type){
        if(cigarS.isEmpty()){
            cigarS.addFirst(new Tuple<>(1,type));
        }

        Tuple<Integer,Character> curr = cigarS.get(0);
        if(curr.item2 == type){
            // když je type stejný, tak jenom zvyšujeme velikost
            cigarS.get(0).item1++;
        }
        else{
            cigarS.addFirst(new Tuple<>(1,type));
        }
    }

    /**
     * Vypíše alignment rozdělený do stejně velkých bloků definovaných lineLen. <br>
     * Využití pro vypisování alignmentu do souboru. <br>
     * Posupně všechny znaky zapisuje do řetězce velikosti lineLen, které
     * jsou pak spárované v Tuple. V item1 je referenční sekvence, v item2
     * je porovnávaná sekvence.
     *
     * @param lineLen počet znaků v jednom řetězci
     * @return pole Tuplů se stringy délky lineLen
     */
    Tuple<String,String>[] toStringArrays(int lineLen){
        int nLines = alignment.size()/lineLen;
        Tuple<String,String>[] lines = new Tuple[nLines];

        StringBuilder stringRef = new StringBuilder();
        StringBuilder stringComp = new StringBuilder();

        for(int i=0; i<nLines;i++){
            for(int j=0;j<lineLen;j++){
                Tuple<Character, Character> t = alignment.get(i*j);
                stringRef.append(t.item1);
                stringComp.append(t.item2);
            }
            lines[i] = new Tuple<>(stringRef.toString(),stringComp.toString());
            // reset stringBuilderů, aby nebylo nutné pokaždé tvořit novou instanci
            stringRef.setLength(0); stringComp.setLength(0);
        }
        return lines;
    }

    /**
     * Konverze cigarS - ve formátu List<<\Tuple<\Integer,Character>> do jednoho Stringu. <br>
     * Jenom se spojí int a char vedle sebe a řetězí se s dalšími v seznamu.
     *
     * @param lineLen délka jedné řádky řetězce
     * @return pole Stringů dělky lineLen
     */
    String[] cigarToStringArray(int lineLen){
        int nLine = cigarS.size()/lineLen;
        String[] lines = new String[nLine];
        StringBuilder stringBuilder = new StringBuilder();

        int i = 0;
        for(Tuple<Integer,Character> t:cigarS){
            stringBuilder.append(t.item1);
            stringBuilder.append(t.item2);
            if(stringBuilder.length() >= lineLen){
                lines[i++] = stringBuilder.toString();
                stringBuilder.setLength(0);
            }
        }
        return lines;
    }
}
