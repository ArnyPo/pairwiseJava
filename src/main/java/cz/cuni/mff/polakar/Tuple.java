package cz.cuni.mff.polakar;

public class Tuple<T1,T2> {
    T1 item1;
    T2 item2;

    /** Třída pro dvojici definovaá na T1 a T2
     *
     * @param item1 Objekt 1
     * @param item2 Objekt 2
     */
    public Tuple(T1 item1, T2 item2){
        this.item1 = item1;
        this.item2 = item2;
    }
    T1 getItem1(){
        return item1;
    }
    T2 getItem2(){
        return item2;
    }
    void update(T1 item1, T2 item2){
        this.item1 = item1;
        this.item2 = item2;
    }
    @Override
    public String toString(){
        return "(" + item1 + "," + item2 + ")";
    }

    public boolean equals(T1 t1,T2 t2){
        return item1.equals(t1) && item2.equals(t2);
    }
}
