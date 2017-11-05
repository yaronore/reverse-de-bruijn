package reverse;
import gurobi.*;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;


// Represents an edge in the graph.
class Edge
{
	public Vertex     source;
    public Vertex     dest;   // Second vertex in Edge
    int flow, capacity, cost, reducedCost;
    
    public Edge(Vertex s, Vertex d)
    {
    	source = s;
        dest = d;
        flow = capacity = 0;
        cost = 1;
    }
    
    public String getName() {
    	return source.name + dest.name.substring(dest.name.length()-1);
    }
}


// Represents a vertex in the graph.
class Vertex implements Comparable<Vertex>
{
    public String     name;   // Vertex name
    public List<Edge> adj;    // Adjacent vertices
    public boolean    finished; 
    public int distance;

    public Vertex( String nm )
      { name = nm; adj = new LinkedList<Edge>(); finished = false;}
    
    public int getPaliEdgo() {
    	int size = adj.size();
    	for (int i = 0; i < size; i++) {
   		if (isPali(adj.get(i).getName()))
    			return i;
    	}
    	return -1;
    }
    
    public static boolean isPali(String name) {
    	int half = name.length() / 2;
    	if (name.substring(0, half).compareTo(reverse(name.substring(half+(name.length()%2)))) == 0)
    		return true;
    	else
    		return false;
    }
    
    private static String reverse(String name) {
    	String rc="";
    	for (int i = name.length()-1; i >=0; i--)
    		rc += name.charAt(i);
    	return rc;
    }

    @Override public int compareTo(Vertex y) {
           // Assume neither string is null. Real code should
            // probably be more robust
            if (distance < y.distance)
            {
                return -1;
            }
            if (distance > y.distance)
            {
                return 1;
            }
            return 0;
        }
}

public class Graph
{
    private Vertex []vertexArray; // Arrays of vertices
    private int numEdges;		  // Number of edges
    private int k;				  // k
    private int alphabetSize;
    private String alphabet;
    
    public Graph(String _alphabet, int _k, boolean flow) {
    	k = _k;
    	alphabet = _alphabet;
  
    	// Allocate array of vertices
    	int add = 0; if (flow) add = 2;
    	alphabetSize = alphabet.length();
    	vertexArray = new Vertex[(int)Math.pow(alphabetSize, k-1)+add];
    	for (int i = 0; i < vertexArray.length; i++)
    		vertexArray[i] = null;
    }
    
    // Return the number of edges
    public int numEdges() {
    	return numEdges;
    }
    
    // Add edges and initialize vertices
    public void generateGraph() {

    	// Insert all (k-1)-mers as vertices
		for (int i=0; i < Math.pow(alphabetSize, k-1); i++) {
			String istr = getString(i, k-1);
			for (int j = 0; j < alphabetSize; j++) {
				Edge e = addEdge(istr, istr.substring(1) + getChar(j));
			}
		}
    }
    
    private void checkGraph() {
		int []check = new int[(int)Math.pow(alphabetSize,k)];
		for (int i = 0; i < check.length; i++)
			check[i] = 0;
		
		int num = 0;
		for (int i = 0; i < vertexArray.length; i++) {
			if (vertexArray[i].adj.size() < alphabetSize-1)
			System.out.println("Degree = " + vertexArray[i].adj.size() + " " + vertexArray[i].name);
			for (int j = 0; j < vertexArray[i].adj.size(); j++) {
				String source = vertexArray[i].adj.get(j).source.name;
				String dest = vertexArray[i].adj.get(j).dest.name;
				if (source.length() > 1 && dest.length() > 1) {
					if (source.substring(1,k-1).compareTo(dest.substring(0,k-2)) != 0) {
						num++;
					}
					else
						check[getInt(source + dest.substring(k-2))]++ ;
				}
			}
		}
		
		int[] outdegree = new int[vertexArray.length];
		for (int i = 0; i < vertexArray.length; i++)
			for (int j = 0; j < vertexArray[i].adj.size(); j++)
				outdegree[getInt(vertexArray[i].adj.get(j).dest.name)]++;
		for (int i = 0; i < vertexArray.length; i++) {
			if (outdegree[i] != vertexArray[i].adj.size())
				System.out.println("Unbalanced " + i + " " + vertexArray[i].name + " " + outdegree[i] + " " + vertexArray[i].adj.size());
			else if (Vertex.isPali(vertexArray[i].name) && outdegree[i] % 2 == 1)
				System.out.println("Odd degree for palindromic " + i + " " + vertexArray[i].name);
			if (Vertex.isPali(vertexArray[i].name))
				System.out.println(vertexArray[i].name + " " + outdegree[i] + " " + vertexArray[i].adj.size());
		}
		
//		System.out.println("A total of " + num + " don't agree");
		
		for (int i = 0; i < check.length; i++) {
			if (check[i] > 1) System.out.println("More than one edge: " + getString(i, k) + " " + check[i]);
			if (check[i] - check[getInt(getRev(getString(i, k)))] != 0)
				System.out.println("Rev has different number: " + check[i] + " " + check[getInt(getRev(getString(i, k)))] + " " + getString(i, k));
		}
		
    }
    
    // Add a new edge to the graph
    private Edge addEdge( String sourceName, String destName)
    {
        Vertex v = getVertex( sourceName );
        Vertex w = getVertex( destName );
        Edge e = new Edge(v,w);
        v.adj.add(0, e );
        numEdges++;
        return e;
    }
    
    // Get edge
    private Edge getEdge(String sourceName, String destName)
    {
       Vertex v = getVertex( sourceName );
       Vertex w = getVertex( destName );
      	for (int i = 0; i < v.adj.size(); i++)
  		if (v.adj.get(i).dest.name.compareTo(w.name) == 0) {
        			return v.adj.get(i);
        } 
      return null;
    }
    
    private boolean removeEdge(String kmer) 
    {
        Vertex v = getVertex( kmer.substring(0,kmer.length()-1));
        Vertex w = getVertex( kmer.substring(1,kmer.length()));
    	boolean removed = false;
    	for (int i = 0; i < v.adj.size() && !removed; i++)
    		if (v.adj.get(i).dest.name.compareTo(w.name) == 0) {
    			v.adj.remove(i); removed = true; numEdges--;
    		}
    	return removed;
    }
    
    private Edge addEdge(String kmer)
    {
        Vertex v = getVertex( kmer.substring(0,kmer.length()-1));
        Vertex w = getVertex( kmer.substring(1,kmer.length()));
        Edge e = new Edge(v,w);
        v.adj.add(0, e );
        numEdges++;
        return e;
    }

    // Get vertex, if does not exist, allocate it
    private Vertex getVertex( String vertexName )
    {
        Vertex v = vertexArray[getInt(vertexName)];
        if( v == null )
        {
            v = new Vertex( vertexName );
            vertexArray[getInt(vertexName)] = v;
        }
        return v;
    }
      
    // The main algorithm to find Eulerian cycle
	public String findEuler(int k) {

		// Keep all paths in a stack, later attach them
		Stack<LinkedList<Edge> > allPaths = new Stack<LinkedList<Edge> >();
		Stack<Integer> startingPoints = new Stack<Integer>();
		
		// The path worked on
		LinkedList<Edge> path = new LinkedList<Edge>();
		startingPoints.add(0);
		
		// The current vertex and edge
		String startingVertex = "";
		for (int i = 0; i < k - 1; i++) 
			startingVertex += alphabet.charAt(0);
		Vertex curr = getVertex(startingVertex);
		Edge e;
		
		// Random and indices
		int i;
		int prevIndex = -1;
		
		// Continue while not all edges have been traveresed
		while (numEdges > 0) {
			
			// Check if vertex has any edges to traverse
			if (curr.adj.size() == 0) {

				// Vertex has no edges to continue from
				curr.finished = true;

				// Need to start looking for other vertices, from this path or previous
				if (prevIndex == -1) i = path.size()-1;
				else i = prevIndex;

				// Look for unfinished vertex
				do {
					curr = path.get(i).source;
					if (curr.adj.size() == 0) curr.finished = true;
					i--;
				} while (curr.finished && i >= 0);
				
				// Found no unfinished vertex, thus continue from previous path
				if (i == -1 && curr.finished) {
					// First, attach this path to previous path
					LinkedList<Edge> prevPath = allPaths.pop();
					int l = startingPoints.lastElement();
					prevPath.addAll(startingPoints.pop()+1, path);
					path = prevPath;
					curr = path.getLast().dest;
					
					// Remember to continue from this index to search
					prevIndex = l;

				} else {
					// Found unfinished vertex, remember current path
					startingPoints.push(i);
					allPaths.push(path);
					prevIndex = -1;
					path = new LinkedList<Edge>();
				}
			}

			else {
				// Vertex has edges to traverse, pick edge arbitrarily
				i = 0;
				e = curr.adj.get(i);
				path.add(e);
				curr.adj.remove(e);
				numEdges--;
			
				// Remove the reverse edge
				Vertex rcdest = getVertex(getRev(curr.name));
				curr = e.dest;
				Vertex rcsrc = getVertex(getRev(e.dest.name));
				
				e = null;
				i = 0;
				do {
					e = rcsrc.adj.get(i);
					i++;
				}while (e.dest.name.compareTo(rcdest.name)!= 0);
				rcsrc.adj.remove(e);
				numEdges--;
			}
		}
		// Push last path to the stack
		allPaths.push(path);
		
		// Attach path, starting from the last to the first
		int size = allPaths.size();
		
		for (i=0; i < size-1; i++) {
			path = allPaths.pop();
			LinkedList<Edge> prevPath = allPaths.pop();
			int l = startingPoints.pop();
			System.out.println("starting point " + i + " " + l);
			prevPath.addAll(l+1, path);
			allPaths.push(prevPath);
		}
	
//		System.out.println("  Attached all paths.");
		// Create the sequence from the path
		LinkedList<Edge> totalPath = allPaths.get(0);
		char []seq = new char[totalPath.size()];
		i = 0;
		while (totalPath.size() > 0) {
	//		System.out.println("total path " + totalPath.getFirst().source.name);
			seq[i++] = totalPath.getFirst().source.name.charAt(k-2);
			totalPath.removeFirst();
		}
		
//		System.out.println("  Created sequence.");

		return new String(seq);
	}


	

	


	
	public static boolean isHomo(String a, String b) {
		char first = a.charAt(0);
		for (int i = 0; i < a.length(); i++) {
			if (a.charAt(i) != first) return false;
		}
		for (int i = 0; i < b.length(); i++) {
			if (b.charAt(i) != first) return false;
		}
		return true;
	}
	
	public static boolean isHomo(String a) {
		char first = a.charAt(0);
		if (first == 's' || first == 't') return false;
		for (int i = 0; i < a.length(); i++)
			if (a.charAt(i) != first) return false;
		return true;
	}
	
    
    private String getRev(String str) {
    	String result = "";
    	for (int i = str.length()-1; i >= 0; i--) {
    		result += str.charAt(i);
    	}
    	return result;
    }
   
    
    public void augmentGraph(boolean augment) {
    	
    	if (k % 2 == 0) {
    		// Remember which polynomials were handled
    		Hashtable<String, Integer> set = new Hashtable<String, Integer>();
    	
    		// Iterate over all polynomials
    		for (int i = 0; i < Math.pow(alphabetSize, k/2); i++) {
    			String half = getString(i, k/2);		
    			String curr = half + getRev(half);
    			
    			// Add edges
        		int num = 0;
    			for (int l = 0; l < k; l++) {
    				String src = curr.substring(0,k-1);
    				String dest = curr.substring(1,k);
    				
    				// Add only if the edges were not inserted
    				if (!set.containsKey(src+dest)) {
    					if (augment) addEdge(curr); else removeEdge(curr);
    					set.put(src+dest, 0);num++;
    				}
    				curr = curr.substring(1,k) + curr.charAt(0);
    			}

    		}
    	}
    	else {
       		// Remember which polynomials were handled
       		Hashtable<String, Integer> set = new Hashtable<String, Integer>();
        	
       		// Iterate over all polynomials
       		for (int i = 0; i < Math.pow(alphabetSize, (k+1)/2); i++) {
       			String half = getString(i, (k+1)/2);
       			String curr = half + getRev(half).substring(1);
       			
        			// Add edges
            		int num = 0;
        			for (int l = 0; l < k; l++) {
        				String src = curr.substring(0,k-1);
        				String dest = curr.substring(1,k);
        				
        				// Add only if the edges were not inserted
        				if (!set.containsKey(src+dest)) {
        					if (augment) addEdge(curr); else removeEdge(curr);
        					set.put(src+dest, 0);num++;
        				}
        				curr = curr.substring(1,k) + curr.charAt(0);
        			}
        		}
    	}
    	
    	// check that all palindromic vertices have even degrees
    	HashSet<String> set = new HashSet<String>();
    	for (int i = 0; i < vertexArray.length; i++) {
    		if (Vertex.isPali(vertexArray[i].name) && vertexArray[i].adj.size() %2 == 1)
    			set.add(vertexArray[i].name);
    	}
    	String[] vertices = set.toArray(new String[0]);
    	for (int i = 0; i < vertices.length; i+=2) {
			String str1 = vertices[i];
			String str2 = vertices[i+1];
			String str = str1 + str2 + str1;
			
			for (int j = 0; j < 2*(k-1); j++) {
				addEdge(str.substring(j, j+k));
			}
    	}
    }
    
    public boolean checkRev(String seq, int k) {
    	// Create counter for each k-mer
		int []check = new int[(int)Math.pow(alphabetSize,k)];
		for (int i = 0; i < check.length; i++)
			check[i] = 0;
		
		// Count how many times each k-mer or its RC appeared
		for (int i = 0; i < seq.length()-k+1; i++) {
			check[getStrInt(alphabetSize, seq.substring(i, i+k))]++;
			check[getStrInt(alphabetSize, getRev(seq.substring(i,i+k)))]++;
		}
		
		boolean result = true;
		for (int i = 0; i < check.length; i++)
			if (check[i] < 1) result = false;
		
		return result;
	}
	
	private String getString(int i, int k) {
		return getStrString(i, k);
	}
	
	private String getStrString(int i, int k) {
		String result = "";
		while (result.length() < k) {
			result += getChar(i % alphabetSize);
			i = i / alphabetSize;
		}
		return result;
	}
	
	private char getChar(int i) {
		return alphabet.charAt(i);
	}
	
	private int getInt(String str) {
		return getStrInt(alphabetSize, str);
	}
	
    private int getStrInt(int alphabetSize, String str) {
    	int result = 0;
    	
    	for (int i = 0; i < str.length(); i++)
    		result += getOneInt(str.charAt(i)) * Math.pow(alphabetSize, i);
    	
    	return result;
    }
    private int getOneInt(char a) {
    	int i = 0;
    	while (alphabet.charAt(i) != a) i++;
    	return i;
    }
    
    public void ILP(int time, String seq) throws GRBException {
        
        GRBEnv env = new GRBEnv("reverse_"+k+"_"+alphabetSize+".log");
        env.set(GRB.DoubleParam.TimeLimit, time);
        GRBModel model = new GRBModel(env);
        
        // Create variables
        int nume = (int)Math.pow(alphabetSize, k);
        int numv = (int)Math.pow(alphabetSize, k-1);
        GRBVar[] kmers = model.addVars(nume, GRB.INTEGER); // Indicator variables if vertex is removed
        model.update();
       
        // Objective
        GRBLinExpr expr = new GRBLinExpr();
        double[] coeffs = new double[nume];
        Arrays.fill(coeffs, 1.0);
        expr.addTerms(coeffs, kmers);
        model.setObjective(expr, GRB.MINIMIZE);
        
        // Set feasible solution
        int[] count = new int[nume];
        for (int i = 0; i < seq.length()-k+1; i++) {
        	String kmer = seq.substring(i,  i+k);
        	count[getInt(kmer)]++;
        }
        for (int i = 0; i < nume; i++) {
        	expr = new GRBLinExpr();
            kmers[i].set(GRB.DoubleAttr.Start, count[i]);
    	}
    	
        // Coverage constraint
        for (int i = 0; i < nume; i++) {
        	if (i < getInt(getRev(getString(i,k)))) {
        	expr = new GRBLinExpr();
        	expr.addTerm(1.0, kmers[i]);
        	expr.addTerm(1.0, kmers[getInt(getRev(getString(i,k)))]);
        	model.addConstr(expr, GRB.GREATER_EQUAL, 1, "cov"+i);
        	}
        	if (i == getInt(getRev(getString(i,k)))) {
            	expr = new GRBLinExpr();
            	expr.addTerm(1.0, kmers[i]);
            	model.addConstr(expr, GRB.GREATER_EQUAL, 1, "cov"+i);
            	}
        	}
        
        // Sequence constraint
        for (int i = 0; i < numv; i++) {
        	String kmer = getString(i, k-1);
        	expr = new GRBLinExpr();
        	for (int j = 0; j < alphabet.length(); j++) {
        		expr.addTerm(1.0, kmers[getInt(alphabet.charAt(j)+kmer)]);
        		expr.addTerm(-1.0, kmers[getInt(kmer+alphabet.charAt(j))]);
        	}
        	model.addConstr(expr, GRB.EQUAL, 0, "seq"+i);
        }

        model.update();
        model.optimize(); 
        
        double[] xvals = model.get(GRB.DoubleAttr.X, model.getVars());
        for (int i = 0; i < nume; i++) {
        	if (i < getInt(getRev(getString(i,k)))) {
            	for (int j = 0; j < xvals[i]+xvals[getInt(getRev(getString(i,k)))]-1; j++) {
            		addEdge(getString(i,k));
            		addEdge(getRev(getString(i,k)));
            	}
            }
        	if (i == getInt(getRev(getString(i,k)))) {
        		for (int j = 0; j < xvals[i]*2-1; j++) {
            		addEdge(getString(i,k));
            	}
        	}
        }
        
        // Dispose of model and environment
        model.dispose();
        env.dispose();                   
        
        return;
    }
}


