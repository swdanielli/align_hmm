package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TextbookClass {
	//public String truth [];
	public int observationId [][];
    public double observationCount [][];
    public String sectionId [];
    /*
	public int usedId [];
	public int cvSetSeed;
	public int cvFlag;
	public int cvFold;
	public boolean isTimeIndex;*/
	
	public TextbookClass(int observationId[][], double observationCount[][],
			String sectionId[]) {
		this.observationId = observationId;
	    this.observationCount = observationCount;
	    this.sectionId = sectionId;
	}
	
	/** TextbookClass
	 * 
	 * @param textbookFileName
	 *            o[t], sentence or segment based, line split, sparse
 	 */
	public TextbookClass(String textbookFileName) {
		loadTextbook(textbookFileName);
	}
	
	/** load textbook
	 * 
	 * @param textbookFileName
	 *            o[t], sentence or segment based, line split, sparse
 	 */
	public void loadTextbook (String textbookFileName) {
		List<int[]> arrayArrayBuffer = null;
        List<double[]> arrayArrayCntBuffer = null;
        List<String> sectionIdBuffer = new ArrayList<String>();
        		
		BufferedReader br = null;
		/* Load textbook, o */
		try { br = new BufferedReader(new FileReader(textbookFileName)); }
		catch (FileNotFoundException e) { e.printStackTrace(); }
        
		try {
			String line = br.readLine();
			arrayArrayBuffer = new ArrayList<int[]>();
			arrayArrayCntBuffer = new ArrayList<double[]>();
			String tbBuffer = "";
			while (line != null) {
				String tags[] = line.split("\t");
				
				tbBuffer = tags[0].trim();
				int[] arrayBuffer = null;
				double[] arrayCntBuffer = null;

				if (tbBuffer.length() > 0) {
					String wordId[] = tbBuffer.split("\\s+");
					arrayBuffer = new int[wordId.length];
					arrayCntBuffer = new double[wordId.length];

					for (int i = 0; i < wordId.length; i++) {
						String slots[] = wordId[i].split(":");
						arrayBuffer[i] = Integer.parseInt(slots[0]);
						arrayCntBuffer[i] = Float.parseFloat(slots[1]);
					}
				} else {
					arrayBuffer = new int[0];
					arrayCntBuffer = new double[0];
				}

				arrayArrayBuffer.add(arrayBuffer);
				arrayArrayCntBuffer.add(arrayCntBuffer);
				sectionIdBuffer.add(tags[1]);
				line = br.readLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	    
		sectionId = new String[sectionIdBuffer.size()];
		for (int i = 0; i < sectionIdBuffer.size(); i++)
			sectionId[i] = sectionIdBuffer.get(i);
		
		observationId = new int[arrayArrayBuffer.size()][];
		for (int i = 0; i < arrayArrayBuffer.size(); i++)
			observationId[i] = arrayArrayBuffer.get(i);

		observationCount = new double[arrayArrayCntBuffer.size()][];
		for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
			observationCount[i] = arrayArrayCntBuffer.get(i);

	    /*
	    for (int i = 0; i < observation.length; i++) {
			for (int j = 0; j < observation[i].length; j++)
				System.out.print(observation[i][j] + " ");
			System.out.println("");
		}
	    */
	}
}
