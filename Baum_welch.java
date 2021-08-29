package Assignment_2;

import Jama.*; 

public class Baum_welch {
	 int a[][];
	 int b[][];
	public Matrix forward(int[] V,double[][] a1,double[][] b1,double[][] initial_distribution1){
		Matrix a = new Matrix(a1);
		Matrix b = new Matrix(b1);
		double[][] alpha1 = new double[V.length][a.getRowDimension()];
		Matrix alpha = new Matrix(alpha1);
		Matrix initial_distribution = new Matrix(initial_distribution1);
		alpha.setMatrix(0,0,0,1,(initial_distribution.transpose()
				.arrayTimesEquals(b.getMatrix(0,b.getRowDimension()-1,V[0],V[0])).transpose()));
		for (int t=1;t<V.length;t++){
			for(int j=0; j<a.getRowDimension();j++){
				Matrix f = alpha.getMatrix(t-1,t-1,0,alpha.getColumnDimension()-1).transpose()
						.arrayTimesEquals(a.getMatrix(0,a.getRowDimension()-1,j,j));
				alpha.set(t,j,f.norm1()* b.get(j, V[t]));
			}
		}
		return alpha;
	}
	
	public Matrix backward(int[] V,double[][] a1,double[][] b1){
		Matrix a = new Matrix(a1);
		Matrix b = new Matrix(b1);
		int len =V.length;
	    Matrix ones = new Matrix(a.getRowDimension(), 1);
	    for (int i = 0; i < a.getRowDimension(); i++)
	      ones.set(i, 0, 1.0);
		Matrix beta = new Matrix(V.length,a.getRowDimension());
		beta.setMatrix(V.length-1,V.length-1,0,1,ones.transpose());
		
		for (int t=V.length-2;t>=0;t = t-1){
			for(int j=0; j<a.getRowDimension();j++){
				Matrix f = beta.getMatrix(t+1,t+1,0,beta.getColumnDimension()-1).transpose().
						arrayTimes(b.getMatrix(0,b.getRowDimension()-1,V[t+1],V[t+1]));
				Matrix f1 = f.transpose().arrayTimes(a.getMatrix(j,j,0,a.getRowDimension()-1));
		
				double sum = 0.0;
				for (int i = 0; i < 1; i++) {
					for (int k = 0; k < 2; k++) {
					sum += f1.get(i, k);
					}
				}
				beta.set(t,j,sum);
			}
		}
		return beta;
	}
	 double a_[][] = {{3.18, 6.818},{2.3403, 1.000}};
	 Matrix a1 = new Matrix(a_);
	 double b_[][] = {{3.134, 0.000, 1.000},{3.817, 2.545, 3.63}};
	 Matrix b1 = new Matrix(b_);
	public int dotProduct(int vect_A[], int vect_B[])
	    {
	        int product = 0;
	        int n = vect_A.length;
	        // Loop for calculate cot product
	        for (int i = 0; i < n; i++)
	            product = product + vect_A[i] * vect_B[i];
	        return product;
	    }
	 public void baum_welch(int[] V,double[][] a1, double[][] b1, double[][] initial_distribution,int n_iter){
		 int M = a1.length;
		 int T = V.length;
		 Matrix a = new Matrix(a1);
		 Matrix b = new Matrix(b1);
		 for(int n=0;n<n_iter;n++){
			 Matrix alpha= forward(V, a1, b1, initial_distribution);
			 Matrix beta = backward(V, a1, b1);
		 int[][][] xi = new int[M][M][T-1];
		 for(int t=0;t<T-1;t++){
			 Matrix o = alpha.getMatrix(t+1,t+1,0,alpha.getColumnDimension()-1).transpose().arrayTimes(a);
			 int sum = 0;
			 for (int i1 = 0; i1 < 1; i1++) {
					for (int k1 = 0; k1 < 2; k1++) {
					sum += o.get(i1, k1);
					}
				}
			 Matrix o1 = b.getMatrix(0,b.getRowDimension()-1,V[t+1],V[t+1]).times(sum);
			 Matrix o2 = o1.times(beta.getMatrix(t+1,t+1,0,beta.getColumnDimension()-1).transpose());
			 int denominator = 0;
			 for (int i1 = 0; i1 < 1; i1++) {
					for (int k1 = 0; k1 < 2; k1++) {
					denominator += o2.get(i1, k1);
					}
				}
			 for( int i=0;i<M;i++){
			//	 numerator = alpha[t, i] * a[i, :] * b[:, V[t + 1]].T * beta[t + 1, :].T
				 Matrix numerator1 = a.getMatrix(i,i,0,a.getColumnDimension()-1).times(alpha.get(t, i));
				 Matrix numerator2 = numerator1.times( b.getMatrix(0,b.getRowDimension()-1,V[t+1],V[t+1])).transpose();
				 Matrix numerator = numerator2.times(beta.getMatrix(t+1,t+1,0,beta.getColumnDimension()-1).transpose());
			 }
		 }
		 int gamma =0;
		 for (int i1 = 0; i1 < 1; i1++) {
				for (int k1 = 0; k1 < 2; k1++) {
				gamma += xi[i1][ k1][i1];
				}
			}
		 }
	  
	 }
	 
	// Driver method 
    public static void main(String[] args) 
    { 
    	int states = 2;
    	int state_syms =3;
    	double[][] a = {{0.5,0.5},{0.5,0.5}};
    	double[][] initial_distribution = {{0.5,0.5}};
    	double[][] b= {{0.111, 0.333, 0.5556},{0.166,0.333,0.5}};
    	int[] V ={2,2,0,0,1,1,2,0,2};
    	
          Baum_welch w=new Baum_welch();
          
          Matrix b1 = new Matrix(b);
          double[][] alpha1 = new double[V.length][2];
  		  Matrix alpha = new Matrix(alpha1);
          Matrix initial_distribution1 = new Matrix(initial_distribution);
          
          alpha.setMatrix(0,0,0,1,(initial_distribution1.transpose().arrayTimesEquals(b1.getMatrix(0,b1.getRowDimension()-1,V[0],V[0])).transpose()));
          Matrix u = w.forward(V,a,b,initial_distribution);
          System.out.println("Forward Matrix :");
          u.print(1, 7);
          Matrix u1 = w.backward(V,a,b);
          System.out.println("backward Matrix :");
          u1.print(1, 7);
          System.out.println("Transmittion Matrix :");
          w.a1.print(1, 7);
          System.out.println("Emission Matrix :");
          w.b1.print(1, 7);
          
    } 
     
    
    
}
