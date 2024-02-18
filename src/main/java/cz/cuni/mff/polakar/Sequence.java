package cz.cuni.mff.polakar;

import java.util.ArrayList;
import java.util.List;

public class Sequence {
    String sequence;
    List<String> seqList = new ArrayList<>();
    String header;
    int size;
    int useLen;

    /** Třída pro uchovávání informací o sekvenci.
     *
     * @param header
     */
    public Sequence(String header){
        this.header = header;
    }
    public Sequence(String header, String sequence){
        // TODO debug!
        this.header = header;
        this.sequence = sequence;
        doSize();
    }

    /**
     * Privátní metoda pro vložení informací do atributů size a useLen
     */
    private void doSize(){
        size = sequence.length();
        useLen = size+1;
    }

    /**
     * Metoda pro ukočení práce s sekvencí a uložení jednotlivých iformací do atributů
     */
    void end(){
        StringBuilder sb = new StringBuilder();
        for (String s: seqList) {
            sb.append(s.trim());
        }
        sequence = sb.toString();
        doSize();
    }
}
