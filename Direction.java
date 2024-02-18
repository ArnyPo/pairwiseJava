package zapoctovy_program;

public class Direction {
    /**
     *
     */
    int i;
    int j;
    private int score;
    public Direction(int i, int j, int score){
        this.i = i;
        this.j = j;
        this.score = score;
    }
    public Direction(){
        i = 0;
        j = 0;
        score = Integer.MIN_VALUE;
    }
    Tuple<Integer,Integer> getDirTup(){
        return new Tuple<Integer,Integer>(i,j);
    }
    int getScore(){
        return score;
    }
    void modify(int i, int j, int score){
        this.i = i;
        this.j = j;
        this.score = score;
    }

    boolean equals(int i, int j){
        return this.i == i && this.j==j;
    }
}
