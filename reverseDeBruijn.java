package reverse;
import gurobi.GRBException;
import java.io.BufferedWriter;
import java.io.FileWriter;

import gurobi.GRBException;

public class reverseDeBruijn {
	
	public static void main(String []argv) throws GRBException {
		
		// Check if enough arguments were given
		if (argv.length != 4) {
			System.out.println("usage: java -jar reverse.jar <k> <output_file> <alphabet> <time_limit>");
			System.out.println("usage: if time limit is set to 0, no ILP is run");
			// ACGT, ACDEFGHIKLMNPQRSTVWY
			return;
		}
		
		// Parse k from the input
		int k = Integer.parseInt(argv[0]);
		System.out.println("Going to create sequence for " + k + " ");

		String alphabet = argv[2];
		System.out.println("alphabet = " + alphabet);
		// Create a new reverse-de-bruijn graph
		Graph g = new Graph(alphabet, k, false);
		g.generateGraph(); 
		
		/// Augment graph and run Euler tour
		g.augmentGraph(true);
		System.out.println("After augmentation, total number of edges is " + g.numEdges());
		String seq = g.findEuler(k);
		// Add the first k-1 characters to cover all k-mers
		seq += seq.substring(0,k-1);
		
		if (argv[3].compareTo("0") != 0) {
			g = new Graph(alphabet, k, false);
			g.generateGraph(); 
			g.ILP(Integer.parseInt(argv[3]), seq);
			seq = g.findEuler(k);
			// Add the first k-1 characters to cover all k-mers
			seq += seq.substring(0,k-1);
		}
		System.out.println("After augmentation, total number of edges is " + g.numEdges());

		// Write the result to a file
		try{
			FileWriter fstream = new FileWriter(argv[1]);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(seq);
			out.close();
		} catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
		
		// Check the result is indeed reverse-de-bruijn
		boolean semi = g.checkRev(seq, k);

		// Some output
		System.out.println(seq.substring(0,Math.min(seq.length(),20)) + "..." + 
				seq.substring(seq.length()-Math.min(seq.length(),20))+ "\n" +
				"Checked and the sequence is reverse-de-Bruijn: " + semi + "\n" +
				"Its length: " + (seq.length()-(k-1)));
	}
}

