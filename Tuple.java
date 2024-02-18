package zapoctovy_program;

public class Tuple<T1,T2> {
    T1 item1;
    T2 item2;
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
}
