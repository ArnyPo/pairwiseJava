package zapoctovy_program;

import java.util.ArrayList;
import java.util.List;

public class Sequence {
    String sequence;
    List<String> seqList = new ArrayList<>();
    String header;
    int size;
    int useLen;

    public Sequence(String header){
        this.header = header;
    }
    void doSize(){
        size = sequence.length();
        useLen = size+1;
    }
    void end(){
        StringBuilder sb = new StringBuilder();
        for (String s: seqList) {
            sb.append(s.trim());
        }
        sequence = sb.toString();
        doSize();
    }
}
