package beast.phylodynamics.model;

import java.util.ArrayList;
import java.util.Collections;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Description;

/**
 * @author Nicola Felix Mueller 
 * Modified by DAR for pairwise coalescent models
 */
@Description("Saves the Timeseries outside of the Epidemiological class" +
			 "Only used to increase overview in the Epi class")
public class TimeSeries extends CalculationNode {

    public ArrayList<Double> integrationTimes;						// times at which euler integration is performed
    public ArrayList<DoubleMatrix> fSeries; 						// Transmission rates
    public ArrayList<DoubleMatrix> gSeries; 						// Any transition that does not need transmission
    public ArrayList<DoubleMatrix> ySeries;							// population sizes over time 

    public ArrayList<DoubleMatrix> mSeries;							// moving Fractions over time 

    
    public ArrayList<DoubleMatrix> nSTraj;							// Trajectory of susceptibles
    public ArrayList<DoubleMatrix> nITraj;							// Trajectory of infected
    public ArrayList<DoubleMatrix> nRTraj;							// Trajectory of infected
    public ArrayList<DoubleMatrix> nSSTraj;							// Trajectory of SS pairs
    public ArrayList<DoubleMatrix> nSITraj;							// Trajectory of SI pairs
    public ArrayList<DoubleMatrix> nISTraj;							// Trajectory of IS pairs
    public ArrayList<DoubleMatrix> nIITraj;							// Trajectory of II pairs
	
	@Override
	public void initAndValidate() throws Exception {

	}
	
	/**
	 * set Back all the ArrayLists and DoubleMatrices
	 * @param states
	 */
	protected void setBack(int states){
		integrationTimes = new ArrayList<Double>();
		
		fSeries = new ArrayList<DoubleMatrix>();
		gSeries = new ArrayList<DoubleMatrix>();
		ySeries = new ArrayList<DoubleMatrix>();
		mSeries = new ArrayList<DoubleMatrix>();
		nSTraj = new ArrayList<DoubleMatrix>();
		nITraj = new ArrayList<DoubleMatrix>();
		nRTraj = new ArrayList<DoubleMatrix>();
		
		nSSTraj = new ArrayList<DoubleMatrix>();
		nSITraj = new ArrayList<DoubleMatrix>();
		nISTraj = new ArrayList<DoubleMatrix>();
		nIITraj = new ArrayList<DoubleMatrix>();
		
		currS = new DoubleMatrix(states);
		currI = new DoubleMatrix(states);
		currR = new DoubleMatrix(states);
		currN = new DoubleMatrix(states);
		printS = new DoubleMatrix(states);
		printI = new DoubleMatrix(states);
		
		// Added for pairwise models
		currSSPairs = new DoubleMatrix(states,states);
		currSIPairs = new DoubleMatrix(states,states);
		currISPairs = new DoubleMatrix(states,states);
		currIIPairs = new DoubleMatrix(states,states);
		
	}

	
	/*
	 * Add to the ArrayLists
	 */
	
	protected void addS(DoubleMatrix S){
		DoubleMatrix s = new DoubleMatrix();
		s.copy(S);
		nSTraj.add(s);
		s = null;
	}
	
	protected void addI(DoubleMatrix I){
		DoubleMatrix i = new DoubleMatrix();
		i.copy(I);
		nITraj.add(i);
		i = null;
	}
	
	protected void addR(DoubleMatrix R){
		DoubleMatrix r = new DoubleMatrix();
		r.copy(R);
		nRTraj.add(r);
		r = null;
	}
	
	protected void addF(DoubleMatrix F){
		DoubleMatrix f = new DoubleMatrix();
		f.copy(F);
		fSeries.add(f);
		f = null;
	}
	
	protected void addG(DoubleMatrix G){
		DoubleMatrix g = new DoubleMatrix();
		g.copy(G);
		gSeries.add(g);
		g = null;
	}
	
	protected void addY(DoubleMatrix Y){
		DoubleMatrix y = new DoubleMatrix();
		y.copy(Y);
		ySeries.add(y);
		y = null;
	}
	
	protected void addM(DoubleMatrix M){
		DoubleMatrix m = new DoubleMatrix();
		m.copy(M);
		mSeries.add(m);
		m = null;
	}
	
	protected void addSSPairs(DoubleMatrix SS){
		DoubleMatrix ss = new DoubleMatrix();
		ss.copy(SS);
		nSSTraj.add(ss);
		ss = null;
	}
	
	protected void addSIPairs(DoubleMatrix SI){
		DoubleMatrix si = new DoubleMatrix();
		si.copy(SI);
		nSITraj.add(si);
		si = null;
	}
	
	protected void addISPairs(DoubleMatrix IS){
		DoubleMatrix is = new DoubleMatrix();
		is.copy(IS);
		nISTraj.add(is);
		is = null;
	}
	
	protected void addIIPairs(DoubleMatrix II){
		DoubleMatrix ii = new DoubleMatrix();
		ii.copy(II);
		nIITraj.add(ii);
		ii = null;
	}
	

	/*
	 * Get Values From the ArrayLists
	 */
	protected double getT(int i){
		return integrationTimes.get(i);
	}
	
	protected DoubleMatrix getS(int i){
		return nSTraj.get(i);
	}
	
	protected DoubleMatrix getI(int i){
		return nITraj.get(i);
	}
	
	protected DoubleMatrix getR(int i){
		return nRTraj.get(i);
	}
	
	protected DoubleMatrix getF(int i){
		return fSeries.get(i);
	}
	
	protected DoubleMatrix getG(int i){
		return gSeries.get(i);
	}
	
	protected DoubleMatrix getY(int i){
		return ySeries.get(i);
	}
	
	protected DoubleMatrix getM(int i){
		return mSeries.get(i);
	}
	
	protected int size(){
		return nSTraj.size();
	}
	
	
	
	/*
	 * Reverse the ArrayLists
	 */
	protected void reverse(){
		
		Collections.reverse(fSeries);
		Collections.reverse(gSeries);
		Collections.reverse(ySeries);
		Collections.reverse(mSeries);
		Collections.reverse(nSTraj);
		Collections.reverse(nITraj);
		Collections.reverse(nRTraj);
		Collections.reverse(integrationTimes);
		
		Collections.reverse(nSSTraj);
		Collections.reverse(nSITraj);
		Collections.reverse(nISTraj);
		Collections.reverse(nIITraj);
		
	}
	
    public DoubleMatrix currS;
    public DoubleMatrix currI;
    public DoubleMatrix currR;
    public DoubleMatrix currN;
    public DoubleMatrix printI;
    public DoubleMatrix printS;
    
    public DoubleMatrix currSSPairs;
    public DoubleMatrix currSIPairs;
    public DoubleMatrix currISPairs;
    public DoubleMatrix currIIPairs;
	
	/*
	 * Set the Double Matrices
	 */
	protected void setS(DoubleMatrix S){
		currS.copy(S);
	}
	
	protected void setI(DoubleMatrix I){
		currI.copy(I);
	}
	
	protected void setI(int i, double init){
		currI.put(i, init);
	}
	
	protected void setR(DoubleMatrix R){
		currR.copy(R);
	}
	
	protected void setN(DoubleMatrix N){
		currN.copy(N);
	}
	
	protected void setPS(DoubleMatrix S){
		printI.copy(S);
	}
	
	// Added for pairs
	protected void setSSPairs(DoubleMatrix SS){
		currSSPairs.copy(SS);
	}
	
	protected void setSIPairs(DoubleMatrix SI){
		currSIPairs.copy(SI);
	}
	
	protected void setISPairs(DoubleMatrix IS){
		currISPairs.copy(IS);
	}
	
	protected void setIIPairs(DoubleMatrix II){
		currIIPairs.copy(II);
	}
	
	/*
	 * get the Double Matrices
	 */
	protected DoubleMatrix getS(){
		return currS;
	}
	
	protected DoubleMatrix getI(){
		return currI;
	}
	
	protected DoubleMatrix getR(){
		return currR;
	}
	
	protected DoubleMatrix getN(){
		return currN;
	}
	
	protected DoubleMatrix getPS(){
		return printI;
	}
	
	protected DoubleMatrix getPI(){
		return printS;
	}
	
	// Added for pairs
	protected DoubleMatrix getSS(){
		return currSSPairs;
	}
	
	// Added for pairs
	protected DoubleMatrix getSI(){
		return currSIPairs;
	}
	
	// Added for pairs
	protected DoubleMatrix getIS(){
		return currISPairs;
	}
	
	// Added for pairs
	protected DoubleMatrix getII(){
		return currIIPairs;
	}
	


}
