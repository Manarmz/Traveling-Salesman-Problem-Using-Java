# Traveling-Salesman-Problem-Using-Java
This was done by TSP in java group project for instructor: Dr. Garsah Farhan Al-Qarni
## Implementation of TSP in Dynamic Programming
### Dynamic Programming Class Code:
```java
public class dynamic {
    
  private final int N;
  private final int START_NODE;
  private final int FINISHED_STATE;

  private double[][] distance;
  private double minTourCost = Double.POSITIVE_INFINITY;

  private List<Integer> tour = new ArrayList<>();
  private boolean ranSolver = false;

  public dynamic(double[][] distance) {
    this(0, distance);
  }

  public dynamic(int startNode, double[][] distance) {

    this.distance = distance;
    N = distance.length;
    START_NODE = startNode;

    // Validate inputs.
    if (N <= 2) throw new IllegalStateException("TSP on 0, 1 or 2 nodes doesn't make sense.");
    if (N != distance[0].length)
      throw new IllegalArgumentException("Matrix must be square (N x N)");
    if (START_NODE < 0 || START_NODE >= N)
      throw new IllegalArgumentException("Starting node must be: 0 <= startNode < N");
    if (N > 32)
      throw new IllegalArgumentException(
          "Matrix too large! A matrix that size for the DP TSP problem with a time complexity of"
              + "O(n^2*2^n) requires way too much computation for any modern home computer to handle");

    // The finished state is when the finished state mask has all bits are set to
    // one (meaning all the nodes have been visited).
    FINISHED_STATE = (1 << N) - 1;
  }

  // Returns the optimal tour for the traveling salesman problem.
  public List<Integer> getTour() {
    if (!ranSolver) solve();
    return tour;
  }

  // Returns the minimal tour cost.
  public double getTourCost() {
    if (!ranSolver) solve();
    return minTourCost;
  }

  public void solve() {

    // solver
    int state = 1 << START_NODE;
    Double[][] memo = new Double[N][1 << N];
    Integer[][] prev = new Integer[N][1 << N];
    minTourCost = tsp(START_NODE, state, memo, prev);

    // Regenerate path
    int index = START_NODE;
    while (true) {
      tour.add(index);
      Integer nextIndex = prev[index][state];
      if (nextIndex == null) break;
      int nextState = state | (1 << nextIndex);
      state = nextState;
      index = nextIndex;
    }
    tour.add(START_NODE);
    ranSolver = true;
  }

  private double tsp(int i, int state, Double[][] memo, Integer[][] prev) {

    // Done this tour. Return cost of going back to start node.
    if (state == FINISHED_STATE) return distance[i][START_NODE];

    // Return cached answer if already computed.
    if (memo[i][state] != null) return memo[i][state];

    double minCost = Double.POSITIVE_INFINITY;
    int index = -1;
    for (int next = 0; next < N; next++) {

      // Skip if the next city has already been visited.
      if ((state & (1 << next)) != 0) continue;

      int nextState = state | (1 << next);
      double newCost = distance[i][next] + tsp(next, nextState, memo, prev);
      if (newCost < minCost) {
        minCost = newCost;
        index = next;
      }
    }

    prev[i][state] = index;
    return memo[i][state] = minCost;
  }
}
```

## Implementation of TSP in Greedy Algorithm
### Greedy Algorithm Class Code:
```java
public class greedy {
        static double findMinRoute(double[][] n)
    {
        int sum = 0;
        int counter = 0;
        int j = 0, i = 0;
        double min = Integer.MAX_VALUE;
        List<Integer> visitedRouteList
                = new ArrayList<>();

        // Starting from the 0th indexed
        // city i.e., the first city
        visitedRouteList.add(0);
        int[] route = new int[n.length];

        // Traverse the adjacency
        // matrix tsp[][]
        while (i < n.length
                && j < n[i].length) {

            // Corner of the Matrix
            if (counter >= n[i].length - 1) {
                break;
            }

            // If this path is unvisited then
            // and if the cost is less then
            // update the cost
            if (j != i
                    && !(visitedRouteList.contains(j))) {
                if (n[i][j] < min) {
                    min = n[i][j];
                    route[counter] = j + 1;
                }
            }
            j++;

            // Check all paths from the
            // ith indexed city
            if (j == n[i].length) {
                sum += min;
                min = Integer.MAX_VALUE;
                visitedRouteList.add(route[counter] - 1);
                j = 0;
                i = route[counter] - 1;
                counter++;
            }
        }

        // Update the ending city in array
        // from city which was last visited
        i = route[counter - 1] - 1;

        for (j = 0; j < n.length; j++) {

            if ((i != j) && n[i][j] < min) {
                min = n[i][j];
                route[counter] = j + 1;
            }
        }
        sum += min;

        // Started from the node where
        // we finished as well.
        System.out.print("Minimum Cost for greedy : ");
        System.out.println(sum);
        
        return sum;
    }
}

```

## Implementation of TSP in Branch and Bound
### Branch and Bound Class Code:
```java
public class branch {

	//NOTE: N changes based on input size
	static int N =5;

	// final_path[] stores the final solution i.e., the
	// path of the salesman.
	static double final_path[] = new double[N + 1];

	// visited[] keeps track of the already visited nodes
	// in a particular path
	static boolean visited[] = new boolean[N];

	// Stores the final minimum weight of shortest tour.
	static double final_res = Integer.MAX_VALUE;

	// Function to copy temporary solution to
	// the final solution
	static void copyToFinal(int curr_path[])
	{
		for (int i = 0; i < N; i++)
			final_path[i] = curr_path[i];
		final_path[N] = curr_path[0];
	}

	// Function to find the minimum edge cost
	// having an end at the vertex i
	static double firstMin(double adj[][], int i)
	{
		double min = Integer.MAX_VALUE;
		for (int k = 0; k < N; k++)
			if (adj[i][k] < min && i != k)
				min = adj[i][k];
		return min;
	}

	// function to find the second minimum edge cost
	// having an end at the vertex i
	static double secondMin(double adj[][], int i)
	{
		double first = Integer.MAX_VALUE, second = Integer.MAX_VALUE;
		for (int j=0; j<N; j++)
		{
			if (i == j)
				continue;

			if (adj[i][j] <= first)
			{
				second = first;
				first = adj[i][j];
			}
			else if (adj[i][j] <= second &&
					adj[i][j] != first)
				second = adj[i][j];
		}
		return second;
	}

	// function that takes as arguments:
	// curr_bound -> lower bound of the root node
	// curr_weight-> stores the weight of the path so far
	// level-> current level while moving in the search
	//		 space tree
	// curr_path[] -> where the solution is being stored which
	//			 would later be copied to final_path[]
	static void TSPRec(double adj[][], int curr_bound, int curr_weight,
				int level, int curr_path[])
	{
		// base case is when we have reached level N which
		// means we have covered all the nodes once
		if (level == N)
		{
			// check if there is an edge from last vertex in
			// path back to the first vertex
			if (adj[curr_path[level - 1]][curr_path[0]] != 0)
			{
				// curr_res has the total weight of the
				// solution we got
				double curr_res = curr_weight +
						adj[curr_path[level-1]][curr_path[0]];
	
				// Update final result and final path if
				// current result is better.
				if (curr_res < final_res)
				{
					copyToFinal(curr_path);
					final_res = curr_res;
				}
			}
			return;
		}

		// for any other level iterate for all vertices to
		// build the search space tree recursively
		for (int i = 0; i < N; i++)
		{
			// Consider next vertex if it is not same (diagonal
			// entry in adjacency matrix and not visited
			// already)
			if (adj[curr_path[level-1]][i] != 0 &&
					visited[i] == false)
			{
				int temp = curr_bound;
				curr_weight += adj[curr_path[level - 1]][i];

				// different computation of curr_bound for
				// level 2 from the other levels
				if (level==1)
				curr_bound -= ((firstMin(adj, curr_path[level - 1]) +
								firstMin(adj, i))/2);
				else
				curr_bound -= ((secondMin(adj, curr_path[level - 1]) +
								firstMin(adj, i))/2);

				// curr_bound + curr_weight is the actual lower bound
				// for the node that we have arrived on
				// If current lower bound < final_res, we need to explore
				// the node further
				if (curr_bound + curr_weight < final_res)
				{
					curr_path[level] = i;
					visited[i] = true;

					// call TSPRec for the next level
					TSPRec(adj, curr_bound, curr_weight, level + 1,
						curr_path);
				}

				// Else we have to prune the node by resetting
				// all changes to curr_weight and curr_bound
				curr_weight -= adj[curr_path[level-1]][i];
				curr_bound = temp;

				// Also reset the visited array
				Arrays.fill(visited,false);
				for (int j = 0; j <= level - 1; j++)
					visited[curr_path[j]] = true;
			}
		}
	}

	// This function sets up final_path[]
	static void TSP(double adj[][])
	{
		int curr_path[] = new int[N + 1];

		// Calculate initial lower bound for the root node
		// using the formula 1/2 * (sum of first min +
		// second min) for all edges.
		// Also initialize the curr_path and visited array
		int curr_bound = 0;
		Arrays.fill(curr_path, -1);
		Arrays.fill(visited, false);

		// Compute initial bound
		for (int i = 0; i < N; i++)
			curr_bound += (firstMin(adj, i) +
						secondMin(adj, i));

		// Rounding off the lower bound to an integer
		curr_bound = (curr_bound==1)? curr_bound/2 + 1 :
									curr_bound/2;

		// We start at vertex 1 so the first vertex
		// in curr_path[] is 0
		visited[0] = true;
		curr_path[0] = 0;

		// Call to TSPRec for curr_weight equal to
		// 0 and level 1
		TSPRec(adj, curr_bound, 0, 1, curr_path);
	}

}

```
## Main Class With a Five City Problem Array
```java
 public static void main(String[] args) {
        // TODO code application logic here
        long start1 = System.currentTimeMillis();
        double[][] distanceMatrix ={  {0.0,  3.0,  4.0,  2.0,  7.0},
{3.0,  0.0,  4.0,  6.0,  3.0},
{4.0,  4.0,  0.0,  5.0 , 8.0},
{2.0,  6.0,  5.0,  0.0,  6.0},
{7.0,  3.0,  8.0,  6.0,  0.0}    
        }; 
        
            for (int i=0; i<20; i++){
      dynamic D = new dynamic(distanceMatrix);
       D.getTourCost();
        System.out.println("Cost for Dynammic programming: "+ D.getTourCost());
        }
        
        long end1 = System.currentTimeMillis();
       
      float time1 = (end1 - start1);
       System.out.println("The time for dynammic programming :"+ (time1)/20);
       
       
       long start2 = System.currentTimeMillis(); 
       
        double[][] distanceMatrix2 ={ {0.0,  3.0,  4.0,  2.0,  7.0},
{3.0,  0.0,  4.0,  6.0,  3.0},
{4.0,  4.0,  0.0,  5.0 , 8.0},
{2.0,  6.0,  5.0,  0.0,  6.0},
{7.0,  3.0,  8.0,  6.0,  0.0}    
        }; 
        
   for (int i=0; i<20; i++){
        branch b = new branch();
        b.TSP(distanceMatrix);
        System.out.println("Cost for branch and bound: "+b.final_res);
        }
        long end2 = System.currentTimeMillis();
       
      float time2 = (end2 - start2);
       System.out.println("The time for branch and bound : "+ (time2)/20);
      
       
       
              long start3 = System.currentTimeMillis(); 
       
        double[][] distanceMatrix3 ={ {0.0,  3.0,  4.0,  2.0,  7.0},
{3.0,  0.0,  4.0,  6.0,  3.0},
{4.0,  4.0,  0.0,  5.0 , 8.0},
{2.0,  6.0,  5.0,  0.0,  6.0},
{7.0,  3.0,  8.0,  6.0,  0.0}    
        }; 
        
   for (int i=0; i<20; i++){
      greedy G = new greedy();
G.findMinRoute(distanceMatrix3);
        }
        long end3 = System.currentTimeMillis();
       
      float time3 = (end3 - start3);
       System.out.println("The time for greedy : "+ (time3)/20);
  
    }
```
#Appendix
##The Process of Converting TSPLIB Data Into an Array
### Data Reader Class
```java
public static Vertex[] read(File file) {
        Vertex[] arrayVertex = null;
        try {
            FileReader reader = new FileReader(file);
            BufferedReader br = new BufferedReader(reader);
            String line;
            ArrayList<Vertex> listVertex = new ArrayList<>();

            //read header data set
            boolean readVertex = false;
            while ((line = br.readLine()) != null) {
                          
                if (line.equals("NODE_COORD_SECTION")) {
                    readVertex = true;
                } else if (line.equals("EOF")) {
                    readVertex = false;
                    break;
                }

                if (readVertex == true) {
                    line = line.replaceAll("\\s+", " ");
                    if (line.substring(0, 1).equals(" ")) {
                        line = line.substring(1);
                    }
                    String[] dline = line.split("\\s");
                    if (dline.length == 3) {
                        String label = dline[0];
                        double x = Double.parseDouble(dline[1]);
                        double y = Double.parseDouble(dline[2]);
                        Vertex v = new Vertex(label, x, y);
                        listVertex.add(v);
                    }
                }
            }

            if (listVertex.size() > 0) {
                //convert to arrayVertex
                int n = listVertex.size();
                arrayVertex = new Vertex[n];
                for (int i = 0; i < arrayVertex.length; i++) {
                    arrayVertex[i] = listVertex.get(i);
                }
            }

        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return arrayVertex;
    }
}
```
### Vertex Class
```java
public class Vertex {
       
    String label;
    double x;
    double y;

    public Vertex(String label, double x, double y) {
        this.label = label;
        this.x = x;
        this.y = y;
    }

    @Override
    public String toString() {
        return label + "(" + x + "," + y + ")";
    } 
}
```
### Data Class (Calculate Distance)
```java
public class Data {
private Vertex[] arrayVertex = null;
    
public Data (Vertex[] arrayVertex) {
    this.arrayVertex = arrayVertex; 
}

public Vertex[] getArrayVertex () {
    return arrayVertex;
}

public void setArrayvertex (Vertex [] arrayVertex) {
    this.arrayVertex = arrayVertex;
}

    public double calculateDistance(int indexVertex1, int indexVertex2) {
        double distance = -1;
        if (arrayVertex != null
                && indexVertex1 >= 0
                && indexVertex1 < arrayVertex.length
                && indexVertex2 >= 0
                && indexVertex2 < arrayVertex.length) {

            Vertex v1 = arrayVertex[indexVertex1];
            Vertex v2 = arrayVertex[indexVertex2];

            double x1 = v1.x;
            double y1 = v1.y;

            double x2 = v2.x;
            double y2 = v2.y;

            double x = x1 - x2;//selisih x
            double y = y1 - y2;//selisih y

            double xx = x * x;//kuadrat x
            double yy = y * y;//kuadrat y

            double xxyy = xx + yy;
            double z = Math.sqrt(xxyy);
            distance = z;
        }
        return distance;
    }
    
@Override
    public String toString () {
    String result = null;
    if (arrayVertex != null) {
        StringBuffer sb = new StringBuffer();
        sb.append ("Data-------------------------------------------");
        sb.append ("label ( x , y )\n");
        for(int i = 0; i < arrayVertex.length; i++) {
            sb.append(arrayVertex[i].toString()+"\n");
        }
        sb.append (" –----------------------------------------------");
        result = sb.toString () ;
    }
    return result;
}
}

```
### Main Class 
```java
    public static void main(String[] args) {
        // TODO code application logic here
        int x;
        //NOTE: x must change based on input size
        for (x=0; x<18; x++)
        {
        int y=0;
                File fileDataset = new File("a280.tsp");
        Vertex[] arrayVertex = DataReader.read(fileDataset);
        Data myData = new Data(arrayVertex);
        

        
        System.out.print("{ ");
        for(int i=0; i<arrayVertex.length; i++){
            double d = myData.calculateDistance(x, y);
            y++;
       if (i==arrayVertex.length-1 )
             System.out.print(" " + d);
       else
            System.out.print(" " + d+",");
       
    }
        System.out.print("}, \n");
        }
     
}
}
```
