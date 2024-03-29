package cz.cuni.mff.polakar;

public class Match {
    int score;
    Direction dir;
    int gapLen; // TODO používám to vůbec????

    public Match(int score){
        this.score = score;
        this.gapLen = 0;
        dir = new Direction(0,0,score);
    }
    public Match(int score, int i, int j){
        this.score = score;
        gapLen = 0;
        dir = new Direction(i,j,score);
    }
    public Match(Direction dir){
        this.dir = dir;
        score = dir.getScore();
        gapLen = 0;
    }
    public Match(Direction dir, int gapLen){
        this.dir = dir;
        score = dir.getScore();
        this.gapLen = gapLen;
    }
    @Override
    public String toString(){
        return score + " " + dir.toString() + " " + gapLen;
    }
}
