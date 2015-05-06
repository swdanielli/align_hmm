package myHMMV4;

public class tuneConfiguration {
	public int nConfigSet;
	public double bestConfiguration [][];
	public double bestResult [];
	
	public tuneConfiguration(int nConfigSet, double [] bestResult, double bestConfiguration [][]) {
		this.nConfigSet = nConfigSet;
		this.bestResult = bestResult;
		this.bestConfiguration = bestConfiguration;
	}
	
	public tuneConfiguration(int nConfigSet) {
		this.nConfigSet = nConfigSet;
		bestResult = new double [nConfigSet];
		bestConfiguration = new double [nConfigSet][];
		for (int i = 0; i < nConfigSet; i++){
			bestResult[i] = -1.0E10;
			bestConfiguration[i] = null;
		}
	}
	
	public void updateConfiguration(double [] result, double [] configuration) {		
		if (nConfigSet != result.length) System.out.println("Error! Mismatched result length");
		//System.out.println(result[0]);
		for (int i = 0; i < nConfigSet; i++) {
			if (result[i] > bestResult[i]) {
				bestResult[i] = result[i];
				bestConfiguration[i] = configuration;
			}
		}
	}
}
