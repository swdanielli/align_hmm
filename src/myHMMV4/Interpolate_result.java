package myHMMV4;

import java.util.List;

public class Interpolate_result {
	List<TranscriptionClass> transObjs;
    double [] result;
    
    public Interpolate_result(List<TranscriptionClass> transObjs, double [] result) {
            this.transObjs = transObjs;
            this.result = result;
    }
}
