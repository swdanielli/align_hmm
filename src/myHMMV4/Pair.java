package myHMMV4;

public class Pair implements Comparable<Pair> {
	public final int id;
    public final double value;

    public Pair(int id, double value) {
        this.id = id;
        this.value = value;
    }

    @Override
    public int compareTo(Pair other) {
    	/* Descending */
        return -1 * Double.valueOf(this.value).compareTo(other.value);
    }
}
