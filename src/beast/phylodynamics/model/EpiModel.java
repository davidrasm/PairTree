package beast.phylodynamics.model;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;


/**
 * @author David Rasmussen
 */
@Description("Epidemic model for standard random mixing SIR" +
             "Tracks deterministic pop trajectory and F, G and Y matrices")
public class EpiModel extends CalculationNode implements Loggable {
	
    public Input<RealParameter> timeStepInput = new Input<RealParameter>("timeStep",
            "timeStep (defaults to 2).",Input.Validate.REQUIRED);    

    /*
     * Input specifying the origin, meaning the point at which to start with
     * simulating the S(E)IR trajectories as well as the f- g- and y-series 
     */
    public Input<RealParameter> originInput = new Input<RealParameter>("origin",
            "timepoint of origin.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> popNInput = new Input<RealParameter>("popN",
            "nodes in network.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> betaInput = new Input<RealParameter>("beta",
            "transmission rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> nuInput = new Input<RealParameter>("nu",
            "removal rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> R0Input = new Input<RealParameter>("R0",
            "net reproductive rate R0.");
    
    public Input<IntegerParameter> maxKInput = new Input<IntegerParameter>(
    		"maxDegree",
    		"maximum allowed node degree");
    
    /*
     * Boolean Inputs used to define which model to use. Would technically not 
     * be necessary but is nice to have
     */    
    public Input<BooleanParameter> isSeasonalInput = new Input<BooleanParameter>(
    		"isSeasonal",
    		"define weather to use a seasonal forcing model");
    public Input<BooleanParameter> isSIRInput = new Input<BooleanParameter>(
    		"isSIR",
    		"define whether this is an SIS or SIR model");
    public Input<BooleanParameter> isConstantInput = new Input<BooleanParameter>(
    		"isConstant",
    		"define weather to use a constant population model");
    public Input<BooleanParameter> movingFractionInput = new Input<BooleanParameter>(
    		"movingFraction",
    		"define how to calculate the backward migration rate");
    
    public Input<Boolean> estimateR0Input = new Input<Boolean>(
    		"estimateR0",
    		"if estimating R0 and computing tau based on R0. Default is false", false);

    /*
     * Parameters that are needed independent of the population model used
     */
    public boolean dirty;
    public boolean reject = false;
    public boolean diagTrans = false;    
    
    /*
     * TimeSeries Object that stores all the F, G, Y and S I R ArraLists
     */
    public TimeSeries timeS;
    
    /*
     * ArrayList storing all the times between which euler integration is
     * performed
     */
    public ArrayList<Double> integrationTimes;
   
    public int states;
    
    protected double beta;
    protected double nu;
    protected double R0;
    
    protected DoubleMatrix betaMatrix;
   
    /*
     * Boolean that store the population model used for 
     * calculating the F G & Y matrices over time. At least
     * one of those has to be true 
     */
    public boolean seasonal = false;
    public boolean SIR = false;
    public boolean estimateR0 = false;
    public boolean constant = false;
    public boolean movingFraction = false;
    
    private double origin;
    private double timeStep;
    	
	@Override
	public void initAndValidate() throws Exception {
		
		//diagTrans = ratesInput.get().minTransRate;
		//states = ratesInput.get().getStates();
		states = maxKInput.get().getValue();
		timeStep = timeStepInput.get().getValue();
		
		if(originInput.get() != null) origin  = originInput.get().getValue();
		
		/*
		 * Get the population model used as Input
		 */
		if (isSeasonalInput.get() != null) seasonal = isSeasonalInput.get().getValue();
		if (isConstantInput.get() != null) constant = isConstantInput.get().getValue();
		if (isSIRInput.get() != null) SIR = isSIRInput.get().getValue();
		if (movingFractionInput.get() != null) movingFraction = movingFractionInput.get().getValue();
		
		/*
		 * If R0 is estimated than tau will be computed from R0
		 */
		if (estimateR0Input.get() != null) estimateR0 = estimateR0Input.get();

		timeS = new TimeSeries();
	}
	
	public boolean update() throws MathException{
		
		//states = maxKInput.get().getValue() + 1;
		timeS.setBack(states);	// empties all arrays in timeS
		
    	boolean reject = false;
//        if (!dirty) return reject;        
        
        final double endTime = 0.0;
        
        // Pairwise model parameters
        final double popN = popNInput.get().getValue();
        beta = betaInput.get().getValue();
        nu = nuInput.get().getValue();
        R0 = R0Input.get().getValue();
        
        if (estimateR0) {
        	beta = R0 * nu;
        }
        
        // Set up beta matrix
        betaMatrix = DoubleMatrix.zeros(states,states);
        betaMatrix.put(0, 0, beta);
        
        // Set time series at t = 0
        DoubleMatrix initN = new DoubleMatrix(states); initN.put(0,popN);
        DoubleMatrix initI = new DoubleMatrix(states); initI.put(0,1);
        DoubleMatrix initS = initN.sub(initI);
        timeS.setS(initS);
        timeS.setI(initI);
        timeS.setN(initN);
        
        //timeS.setR(timeS.getN().sub(timeS.getS().add(timeS.getI())));
        //timeS.setPS(timeS.getPS());
        //timeS.setPI(timeS.getPI());
        		
        integrationTimes = new ArrayList<Double>();    	
		
		double time = origin;
		
		while (time > endTime){
			integrationTimes.add(time);
			time -= timeStep;
		}		
		integrationTimes.add(endTime);
	
		timeS.addS(timeS.getS());
		timeS.addI(timeS.getI());
		timeS.addR(timeS.getN().sub(timeS.getS().add(timeS.getI())));
		
		int t = 0;
		int dur = integrationTimes.size();		
		do{
			reject = deltaStep(t, dur, timeStep);
			t++;
		} while(integrationTimes.get(t-1)>0 && !reject);

		/*
		 * Reverse the array Lists in order to go from
		 * forward in time to backwards in time
		 */
		timeS.reverse();
		Collections.reverse(integrationTimes);

		return reject;		
	}
	
	/*
	 * Calculation of the exponential trajectories as well as updating of
	 * the initial time point of the first individual
	 */
	protected boolean deltaStep(int t, int dtTimeCount, double TimeStep){
		
		// Transmission
		DoubleMatrix sRowVec = timeS.getS().div(timeS.getN()).transpose();
	    DoubleMatrix transMatrix = betaMatrix.mul(timeS.getI().mmul(sRowVec));
	    DoubleMatrix transVec = transMatrix.columnSums().transpose();  // sum is over rows
		
		// Removal
	    DoubleMatrix recVec = timeS.getI().mul(nu);
	    
	    DoubleMatrix noImmunityVec;
	    if (SIR) {
	    	noImmunityVec = DoubleMatrix.zeros(states);
	    } else {
	    	noImmunityVec = DoubleMatrix.ones(states);
	    }
	    
	    // Change in state variables
	    DoubleMatrix dS = noImmunityVec.mul(recVec).sub(transVec);
	    DoubleMatrix dI = transVec.sub(recVec);
		
		DoubleMatrix Transmission = transMatrix; // transmission is from k --> l
		DoubleMatrix Migration = DoubleMatrix.zeros(states, states); // migration from k --> l
		
		// Have to add F, G and Y
		timeS.addF(Transmission);
		timeS.addG(Migration);
		timeS.addY(timeS.getI());
		
		if (t < (dtTimeCount - 1)) {
		    
	    	// Update state variables
	    	final double dt = integrationTimes.get(t) - integrationTimes.get(t+1);
	    		
	    	//DoubleMatrix newS = timeS.getS().add(dS.mmul(tempDt));
	    	timeS.setS(timeS.getS().add(dS.mul(dt)));
	    	//DoubleMatrix newI = timeS.getI().add(dI.mmul(tempDt));
	    	timeS.setI(timeS.getI().add(dI.mul(dt)));
	    	
	    	// Update pop densities ???
	    	//timeS.setN(timeS.getS().add(timeS.getI()));
	    	   			
	    	if (timeS.getPS().min() < 0|| timeS.getPI().min() < 0){
	    		return true;
	    	}
           	timeS.addS(timeS.getS());
           	timeS.addI(timeS.getI());
           	//timeS.addR(timeS.getR());
           	
           	//DoubleMatrix currI = timeS.getI();
           	//DoubleMatrix currS = timeS.getS();
           	
           	//System.out.println(currI);
           	//System.out.println(currS);

           	
	    }
		
		return false;
	}
	
	
	/*
	private void updateInitialTimePoint(){
		intitalInfected = new DoubleMatrix(states);
		for (int s = 0; s < states; s++)
			intitalInfected.put(s, intitalIntroductionInput.get().getValue(s));
	} */
	
	/*
	 * something else
	 */
	protected boolean delta(int i, int t, double s){
		return false;
	}
	
	public int getTrajLength(){
    	return timeS.size();
    }

    public DoubleMatrix getNumS(int t) {
    	return timeS.getS(t);
    }
    
	public DoubleMatrix getNumI(int t) {
		return timeS.getI(t);
    }
    
	public DoubleMatrix getNumR(int t) {
		return timeS.getR(t);
    }
	
	public double getNumSDeme(int t, int d) {
		return timeS.getS(t).get(d);
    }
    
	public double getNumIDeme(int t, int d) {
		return timeS.getI(t).get(d);
    }
    
	public double getNumRDeme(int t, int d) {
		return timeS.getR(t).get(d);
    }
	
	public double getTime(int t) {
		return integrationTimes.get(t);
    }
    
    /**
     * @param t time interval to return
     * @param k state to return Y for
    */
    public double getY(int t, int k) {
    	if (constant)
    		return timeS.getY(0).get(k);
    	else
    		return timeS.getY(t).get(k);
    }
    
    public DoubleMatrix getY(int t) {
    	if (constant)
    		return timeS.getY(0);
    	else
    		return timeS.getY(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in F
     * @param l column in F
    */
    public double getF(int t, int k, int l) {
    	if (constant)
    		return timeS.getF(0).get(k,l);
    	else
    		return timeS.getF(t).get(k,l);
    }
    
    public DoubleMatrix getF(int t) {
    	if (constant)
    		return timeS.getF(0);
    	else
    		return timeS.getF(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in G
     * @param l column in G
    */
    public double getG(int t, int k, int l) {
    	if (constant)
    		return timeS.getG(0).get(k,l);
    	else
    		return timeS.getG(t).get(k,l);
    }
    
    public DoubleMatrix getG(int t) {
	    if (constant)
			return timeS.getG(0);
		else
			return timeS.getG(t);
    }
    
    /**
     * 
     * @param t
     * @param k
     * @param l
     * @return
     */
    public double getM(int t, int k, int l) {
    	if (constant)
    		return timeS.getM(0).get(k,l);
    	else
    		return timeS.getM(t).get(k,l);
    }
    
    public DoubleMatrix getM(int t) {
	    if (constant)
			return timeS.getM(0);
		else
			return timeS.getM(t);
    }
    
    public Boolean movingFraction(){
    	return movingFraction;
    }
    
    /**
     * @return integrationTimes in reverse order (backwards in time)
     */
    public ArrayList<Double> getIntegrationTimes() {
    	return integrationTimes;
    }
    
    public double getOrigin(){
    	return origin;
    }
    
    public void setOrigin(double origin){
    	this.origin = origin;
    }
    
    /*
     * CalculationNode interface
     */

    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    public void restore() {
        dirty = true;
        super.restore();
    }
    
    public void init(PrintStream out) throws Exception {

        //out.print("meanK\t");
        //out.print("varK\t");
        //out.print("tau\t");

    }

    public void log(int nSample, PrintStream out) {
    	
        //out.format("%g\t", meanK);
    	//out.format("%g\t", varK);
    	//out.format("%g\t", tau);
    	//System.out.println(pk);
        
    	// was this:
        //out.format("%g\t", betaParameter.get().getValue()*n_S_Parameter.get().getValue()
                ///gammaParameter.get().getValue());

        //double tend = NStraj.size() * dt;
        //double delta = tend / (statesToLogInput.get() - 1);

    }

    public void close(PrintStream out) {
    }
    
	
	
}

