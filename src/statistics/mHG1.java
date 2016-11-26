package statistics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class mHG1 {
	private static int threshold=6000;  //The maximal position in the motif binary vector we calculate the p-value for
	private final static int R_ZONE = 0;
	
	
	public static HGScore calculate_HGT(int[] occurrenceVec, int B, int threshold){
		int N = occurrenceVec.length;
		int n = threshold;
		int b = 0; //number of 1's that we have seen so far (i.e. number of 1's from position 0 to position n in the motif_Vector)
		for (int i=0;i<n; i++){
			if (occurrenceVec[i]>0)
				b++;
		}
		return calculate_HGT(N, n, B, b);
	}
	
	public static HGScore calculate_HGT(int N, int[] occurrenceVec, int B){
		int n = occurrenceVec.length;
		int b = 0; //number of 1's that we have seen so far (i.e. number of 1's from position 0 to position n in the motif_Vector)
		for(int occ:occurrenceVec)
			if (occ>0) b++;
		
		return calculate_HGT(N, n, B, b);
	}
	
	public static HGScore calculate_HGT(int N, int n, int B, int b){
		if (B>N || n>N || b>B || b>n) throw new IllegalStateException("Wrong Parameter Sizes!: N="+N+" n="+n+" B="+B+" b="+b);
		double hg = mHG1.HG(N, n, B, b);
		double hgt = mHG1.HGT(hg, n, N, B, b);
		HGScore hgScore = new HGScore(N,B,n,b,hgt);
		return hgScore;
	}
	
	public static HGScore calculate_mHGT(int[] occurrenceVec, int B){
		return calculate_mHGT(occurrenceVec, B, threshold);
	}
	/**
	 * based on Eran's code for mHG. calculates the minimum HG score
	 * by going over all thresholds and returning the optimum.
	 * @param occurrenceVec The ranked list where '1' denotes a member
	 * in the target set.
	 * @param B The size of the target set.
	 * @param threshold The limit of the threshold space.
	 * @return the min HG score.
	 */
	public static HGScore calculate_mHGT(int[] occurrenceVec, int B, int threshold){
		int N = occurrenceVec.length;
		HGScore minHG = new HGScore(N,B,0,0,1);
		if(B==0) return minHG;
		int b = 0; //number of 1's that we have seen so far (i.e. number of 1's from position 0 to position n in the motif_Vector)
		if (B>N) throw new IllegalStateException("Erorr, it is not logical that B is larger than N!: B="+B+" N="+N);
			
		double currHG=1; //HG(N,K,0,0) = 1
		//double mHGT=1;
		long max_index = threshold < N ? threshold : N; 

		for (int n=0; n < max_index; n++){
			if (occurrenceVec[n]>0){
				// get HG(N,B,n+1,b+1) from HG(N,B,n,b)
				currHG = currHG * ((n+1)*(B-b)) / ((N-n)*(b+1));
				b++;
				if (b>B) throw new IllegalStateException("B "+B+" b "+b);
			}else// get HG(N,B,n+1,b) from HG(N,B,n,b)
				currHG = currHG * ((n+1)*(N-B-n+b)) / ((N-n)*(n-b+1));
			
			if (occurrenceVec[n]>0){  //Computing the tail only if we reached 1 in the vector
				if (n+1<max_index && occurrenceVec[n+1]>0) continue; /// the next would be better
				double currHGT = HGT(currHG,(n+1),N,B,b);
				if (currHGT < minHG.score){
					minHG.score = currHGT;
					minHG.b = b;
					minHG.n = (n+1);
				}
			}
			if(b==B) break;
		}

		return minHG;
	}
	
	public static HGScore calculate_mHGT(int[] occurrenceVec, int B, List<Integer> cutoffs, int threshold){
		int N = occurrenceVec.length;
		HGScore minHG = new HGScore(N,B,0,0,1);
		if(B==0) return minHG;
		int b = 0; //number of 1's that we have seen so far (i.e. number of 1's from position 0 to position n in the motif_Vector)
		if (B>N) throw new IllegalStateException("Erorr, it is not logical that B is larger than N!: B="+B+" N="+N);
			
		double currHG=1; //HG(N,K,0,0) = 1
		//double mHGT=1;
		long max_index = threshold < N ? threshold : N; 

		for (int n=0; n < max_index; n++){
			if (occurrenceVec[n]>0){
				// get HG(N,B,n+1,b+1) from HG(N,B,n,b)
				currHG = currHG * ((n+1)*(B-b)) / ((N-n)*(b+1));
				b++;
				if (b>B) throw new IllegalStateException("B "+B+" b "+b);
			}else// get HG(N,B,n+1,b) from HG(N,B,n,b)
				currHG = currHG * ((n+1)*(N-B-n+b)) / ((N-n)*(n-b+1));
			
			/*if (occurrenceVec[n]>0)*/{  //Computing the tail only if we reached 1 in the vector
				if (false == cutoffs.contains(n)) continue;
//				if (n+1<max_index && occurrenceVec[n+1]>0) continue; /// the next would be better
				double currHGT = HGT(currHG,(n+1),N,B,b);
				if (currHGT < minHG.score){
					minHG.score = currHGT;
					minHG.b = b;
					minHG.n = (n+1);
				}
			}
			if(b==B) break;
		}

		return minHG;
	}
	
	public static HGScore calculate_mHGT(int N, int[] occurrenceVec, int B){
		if (B>N) throw new IllegalStateException("Erorr, it is not logical that B is larger than N!: B="+B+" N="+N);

		HGScore minHG = new HGScore(N,B,0,0,1);
		if(B==0) return minHG;
		int b = 0; //number of 1's that we have seen so far (i.e. number of 1's from position 0 to position n in the motif_Vector)
		double currHG=1; //HG(N,B,0,0) = 1
		
		boolean first = true;
		int occ = -1;
		int n=0;
		for (int nextOcc:occurrenceVec){
			if (first){
				occ = nextOcc;
				first = false;
				continue;
			}
			n++;
			if(occ==0){// get HG(N,B,n+1,b) from HG(N,B,n,b)
				currHG = currHG * (n*(N-B-n+b+1)) / ((N-n+1)*(n-b));
				occ = nextOcc;
				continue;
			}
			//if (occ>0)
			b++;
			// get HG(N,B,n+1,b+1) from HG(N,B,n,b)
			currHG = currHG * (n*(B-b+1)) / ((N-n+1)*b);
			
			if (nextOcc==0){
				double currHGT = HGT(currHG,n,N,B,b);
				if (currHGT < minHG.score){
					minHG.score = currHGT;
					minHG.b = b;
					minHG.n = n;
				}
			}
			if(b==B) break;
			
			occ = nextOcc;
		}

		return minHG;
	}
	/**
	 * based on Eran's code. Computes the hyper geometric tail.
	 * @param currHG The current HG(N,K,n,k) score.
	 * @param N The number of elements in vector
	 * @param B The number of 1's in vector
	 * @param n The number of elements so far
	 * @param b The number of 1's in the small vactor that goes until n.
	 * @return The HG tail score
	 */
	public static double HGT(double currHG, int n,int N,int B, int b){
		int min_nB=(n<B) ? n : B;
		double tail=currHG;
		
//		 get HG(N,K,n,i+1) from HG(N,K,n,i)
		for(int i=b;i<min_nB;i++){
			currHG = currHG * ((n-i)*(B-i)) / ((i+1)*(N-n-B+i+1));
			tail += currHG;
		}
		return tail;
	}

	/**
	 * based on Eran's code. Computes the p-value on a mHG score recived
	 * in a B set of N ranked list
	 * @param mHGT The min HG score
	 * @param N The number of elements in vector
	 * @param B The number of 1's in vector
	 * @return The HG tail score
	 */
	public static double calculate_pValue(int N, int B, double mHGT){
		return calculate_pValue(N, B, mHGT, threshold);
	}
	public static double calculate_pValue(int N, int B, double mHGT, int threshold){
		int min_NTHRESHOLD = N < threshold ? N : threshold;
		int min_B_nthresh= B < threshold ? B : threshold; 
		if (B==0) return 1;
		
		double[][] mat = new double[min_B_nthresh+1][min_NTHRESHOLD+1];
		for (int i = 0; i <= min_B_nthresh; i++)
			for (int j = 0; j <= min_NTHRESHOLD; j++)
				mat[i][j] = 0;
		
		mat[0][0] = 1;
		double baseHG = 1;	// holds HG(N,K,n,min(n,K))

		for (int n = 1; n <= min_NTHRESHOLD; n++){
			// n is the number of elemnets in current vector
			int min_nB;
			if (B >= n){
				min_nB = n;
				baseHG = baseHG * (B - n + 1) / (N - n + 1);
			}else{
				min_nB = B;
				baseHG = baseHG * n / (n - B);
			}
			
			if (baseHG <= Double.MIN_VALUE){
				baseHG = Double.MIN_VALUE;
				min_NTHRESHOLD=n;
			}

			double tailHG = baseHG;
			double currHG = baseHG;
			// first loop - sum up the tail, until the sum is bigger than mHGT
			int b;
			for (b = min_nB; tailHG <= mHGT && b > 0; b--){
				// b is the number of ones in current vector
				currHG = currHG * (b*(N-B-n+b)) / ((n-b+1)*(B-b+1));
				//if (currHG == 0) currHG = Double.MIN_VALUE;///
				tailHG += currHG;
				mat[b][n] = R_ZONE;
			}
			// second loop, starts when b is the maximal for which
			// HGT(N,B,n,b)> mHGT
			for (; b > 0; b--){
				// calculate current cell value by two optional cells from
				// which it can be reached
				// 1. last element in vector is 0
				
				mat[b][n]=0; //////////////////////// for priniting reasons
				if (mat[b][n-1]<=1){ //////////////////////// for priniting reasons
					mat[b][n] += mat[b][n-1] * (N-B-n+b+1) / (N-n+1);
				}
				// 2. last element in vector is 1
				if (mat[b-1][n-1]<=1){ //////////////////////// for priniting reasons
					mat[b][n] += mat[b-1][n-1] * (B-b+1) / (N-n+1);
				}
				//if (mat[b][n] == 0) mat[b][n] = Double.MIN_VALUE;///
			}
			mat[0][n] = mat[0][n-1] * (N-B-n+1) / (N-n+1);
			
			if (mat[0][n] == Double.MIN_VALUE){
				min_NTHRESHOLD=n;
				//System.err.println("2: n = "+n);
			}
		}

		double result = 0;
		for (int i = 0; i <= min_B_nthresh; i++)
			result += mat[i][min_NTHRESHOLD];
		return (1 - result);
	}
	
	/**
	 * based on Eran's code. Computes the p-value on a mHG score recived
	 * in a B set of N ranked list. This version of the p-value 
	 * computation uses a box space (spaning W and B) instead of a 
	 * parraralog space (spaning N and B).
	 * @param mHGT The min HG score
	 * @param N The number of elements in vector
	 * @param B The number of 1's in vector
	 * @return The HG tail score
	 */
	public static double calculate_pValue_v2(int N, int B, double mHGT){
		return calculate_pValue_v2(N, B, mHGT, threshold);
	}
	public static double calculate_pValue_v2(int N, int B, double mHGT, int threshold){
		int min_N_nthresh = N < threshold ? N : threshold;
		int min_B_nthresh= B < threshold ? B : threshold;
		//int min_B_nthresh= B;
		
		double[][] mat = new double[min_B_nthresh+1][min_N_nthresh+1];
		for (int i = 0; i <= min_B_nthresh; i++)
			for (int j = 0; j <= min_N_nthresh; j++)
				mat[i][j] = 0;
		
		mat[0][0]=1;
		int W=N-B;
		
		
		double baseHG = 1;	// holds HG(N,K,n,min(n,K))
		double result=0;


		for (int n = 1; n <= min_N_nthresh; n++){
			//1. Calculating the base
			// n is the number of elemnets in current vector
			int min_nB;
			if (n <= B){
				min_nB = n;
				baseHG = baseHG * (B - n + 1) / (N - n + 1);
			}else{
				min_nB = B;
				baseHG = baseHG * n / (n - B);
			}
	 
			double tailHG = baseHG;
			double currHG = baseHG;
			
			//2. From the base going on the diagonal of all the vectors of size n.
			// for each diagonal:  sum up the tail, until the sum is bigger than mHGT
			int b,w;
			for (b = min_nB, w=n-min_nB; tailHG <= mHGT && b >= 0; b--, w++){   ////////////!!!!!!!!! I made b>=0 instead of b>0
				// b is the number of ones in current vector
				currHG = currHG * (b*(N-B-n+b)) / ((n-b+1)*(B-b+1));
				tailHG += currHG;
				mat[b][w] = R_ZONE;
			}
			//Conservation: once this loop is finished the ppoint b w is out of the R_zone space

			//3. going on the remaining part of the diagonal and computing the p_val matirx
			// second loop, starts when b is the maximal for which
			// HGT(N,K,n,k)> mHGT
			if (n==999)
				System.err.println("svcdsfsd");
			for (; b >= 0; b--, w++){
				if ((w==0) && (b==0))
					mat[b][w]=1;
				
				else if ((w==0) && (b>0)) 
					mat[b][w]=mat[b-1][w]*(B-b+1)/(N-n+1);
				
				else if ((w>0) && (b==0)) 
					mat[b][w]=mat[b][w-1]*(W-w+1)/(N-n+1);
				
				else
					mat[b][w]=mat[b-1][w]*(B-b+1)/(N-n+1) + mat[b][w-1]*(W-w+1)/(N-n+1);	
				
	 			if (b+w== min_N_nthresh){ //summing the elements on the diagonal n
					result +=mat[b][w];
				}
			}
		}
		
		return 1-result;
	}

	/**
	 * based on Eran's code. Computes N choose k. The amount of possibilities
	 * to choose 'k' elements from a set of size 'n'.
	 * @param n The universe set size.
	 * @param k The selected set size.
	 * 
	 * @return The amount of possibilities to choose 'k' elements 
	 * from a set of size 'n'.
	 */
	public static long NchooseK(int n, int k)
	{
		if (k < 0 || k > n)
			throw new IllegalArgumentException("n "+n+" k "+k);

		if (k > n - k)
			k = n - k;

		long sum = 1;
		
		int n_k = n - k;
		for(;n > n_k; n--)
			sum *= n;
		for(;k > 1; k--)
			sum /= k;

		return sum;
	}
	
	//	Computes the HG. In order to achieve maximum precision it uses a recursion formula.
	public static double HG(int N, int n_, int B, int b_){
		double currHG=1; //HG for n=0 and b=0
		int n,b;

		if (n_ < b_ || B < b_ || n_ < 0 || b_ < 0){
			System.err.println("Error using HG. n_ "+n_+" must be larger than b_ "+b_+" !!!");
			return -1;
		}
		
		if (n_==0)  //case where n=0, b=0
			return 1;
		double log = 0;
		for (b=0, n=0; b < b_ ; b++, n++)
		{
			double m = (double)((n+1)*(B-b))/((N-n)*(b+1));
			log += Math.log(m);
		}
		
		for (; n < n_; n++)
		{
			double m = (double)((n+1)*(N-B-n+b))/((N-n)*(n-b+1));
			log += Math.log(m);
		}
		
		return Math.exp(log);
		
//		for (b=0, n=0; b < b_ ; b++, n++)
//		{
//			System.out.println(currHG + " b = " + b + " n = " + n);
//			currHG=currHG * ((n+1)*(B-b)) / ((N-n)*(b+1));
//		}
//
//		for (; n < n_; n++)
//			currHG=currHG * ((n+1) * (N-B-n+b)) / ((N-n) * (n-b+1));
		
//		return currHG;
	}
	
	/////////////////////////////////////////////////////////////
//	public static void main(String[] args)throws IOException{
//		String inFileName = "in.txt";
//		String outFileName = "out.txt";
//		BufferedReader in = new BufferedReader(new FileReader(inFileName));
//		List<Integer> binArray = new ArrayList<Integer>();
//		String str;
//		while ((str = in.readLine()) != null) {
//			System.err.println(str);
//			Integer o = Integer.parseInt(str);
//			if (o!=0 && o!=1)
//				throw new IllegalStateException("o "+o+" has to be 0 or 1.");
//			binArray.add(o);
//			
//		}
//		in.close();
//		int[] occVec = new int[binArray.size()];
//		int B = 0;
//		int i=0;
//		for (Integer o:binArray){
//			occVec[i++] = o;
//			if (o == 1)
//				B++;
//		}
//		
//		int threshold = 1000;
//		HGScore score = mHG.calculate_mHGT(occVec, B, threshold);
//		score.calcPvalue(threshold);
//		
//		PrintWriter out = new PrintWriter(outFileName);
//		out.println(score);
//		out.close();
//	}
//}
	
	public static void main(String[] args)throws IOException{
//		int[] a = {0,0,0,0,0,1,1,1,1,0};
//		HGScore s = mHG.calculate_mHGT(10, a, 4);//calculate_HGT(10, 9, 4, 4);
//		System.out.println(s);
//		System.out.println(s.calcPvalue(10));
//		int[] a = {1,1,0,0,0,0,0,0,0,0,0};
//		HGScore s = mHG.calculate_mHGT(11, a, 2);
//		System.out.println(s.calcPvalue(11));
//		int[] b = {1,0,0,0,0,0,0,0,0,0,1};
//		s = mHG.calculate_mHGT(11, b, 2);
//		System.out.println(s.calcPvalue(11));
//		int[] c = {0,0,1,0,0,0,0,0,0,0,1};
//		s = mHG.calculate_mHGT(11, c, 2);
//		System.out.println(s.calcPvalue(11));
//		
//		int[] d = {1,1,0,0,0,0,0,0,0,0,1};
//		s = mHG.calculate_mHGT(11, d, 3);
//		System.out.println(s.calcPvalue(11));
//		
//		int[] e = {1,1,1,0,0,0,0,0,0,0,1};
//		s = mHG.calculate_mHGT(11, e, 4);
//		System.out.println(s.calcPvalue(11));
		
//		HGScore s = mHG.calculate_HGT(4997, 336, 12, 11);
//		System.out.println(s);

	}
}