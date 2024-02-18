package cz.cuni.mff.polakar;

import org.junit.jupiter.api.Test;

class MatrixReaderTest {
    @Test
    void matchMismatchConstructorMatch(){
        MatrixReader mr = new MatrixReader(1,-1,"n");
        assert mr.reader('A','A')==1;
    }
    @Test
    void matchMismatchConstructorMismatch(){
        MatrixReader mr = new MatrixReader(1,-1,"n");
        assert mr.reader('A','C')==-1 ;
    }
}