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
@Description("Pairwise epidemic model simulation for contact networks" +
             "Tracks deterministic pop trajectory and F, G and Y matrices")
public class PairwiseEpiModel extends CalculationNode implements Loggable {
	
    public Input<RealParameter> timeStepInput = new Input<RealParameter>("timeStep",
            "timeStep (defaults to 2).",Input.Validate.REQUIRED);    

    /*
     * Input specifying the origin time of the epidemic
     */
    public Input<RealParameter> originInput = new Input<RealParameter>("origin",
            "timepoint of origin.",Input.Validate.REQUIRED);
    
    /*
     * Input specifying the number of nodes (N) in the network
     */
    public Input<RealParameter> nodesInput = new Input<RealParameter>("nodes",
            "nodes in network.",Input.Validate.REQUIRED);
    
    /*
     * Input specifying the transmission rate
     */
    public Input<RealParameter> tauInput = new Input<RealParameter>("tau",
            "transmission rate.",Input.Validate.REQUIRED);
    
    /*
     * Input specifying the removal (i.e recovery) rate
     */
    public Input<RealParameter> nuInput = new Input<RealParameter>("nu",
            "removal rate.",Input.Validate.REQUIRED);
    
    /*
     * Input specifying the net reproductive number R0 (not used if not est. R0)
     */
    public Input<RealParameter> R0Input = new Input<RealParameter>("R0",
            "net reproductive rate R0.");
    
    /* These now enter only through the degree distribution
    public Input<RealParameter> meanKInput = new Input<RealParameter>("meanK",
            "mean of degree dist.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> varKInput = new Input<RealParameter>("varK",
            "var of degree dist.",Input.Validate.REQUIRED); */
    
    /*
     * Input specifying the maximum degree k in the network
     */
    public Input<IntegerParameter> maxKInput = new Input<IntegerParameter>(
    		"maxDegree",
    		"maximum allowed node degree");
    
    /*
     * Input specifying the degree distribution of the network
     */
    public Input<ParametricDistribution> degreeDistInput = new Input<ParametricDistribution>(
    		"degreeDist",
    		"parametric distribution for degree dist");
    
    /*
     * Boolean Inputs used to define which model to use.
     */    
    public Input<BooleanParameter> isSIRInput = new Input<BooleanParameter>(
    		"isSIR",
    		"define whether this is an SIS or SIR model");
    
    public Input<Boolean> estimateR0Input = new Input<Boolean>(
    		"estimateR0",
    		"if estimating R0 and computing tau based on R0. Default is false", false);

    /*
     * Parameters that are needed independent of the population model used
     */
    public boolean dirty;
    public boolean reject = false;    
    
    /*
     * TimeSeries Object that stores all the F, G, Y as well as state variable trajectories
     */
    public TimeSeries timeS;
    
    /*
     * ArrayList storing integration times (for Euler integration)
     */
    public ArrayList<Double> integrationTimes;
   
    /*
     * Local variables
     */
    public int states;    
    protected double tau;
    protected double nu;
    protected double R0;
 
    public DoubleMatrix pk;
    public DoubleMatrix kSeries;
    public DoubleMatrix qSeries;
    
    public double meanK;
    public double varK;
    
    public DoubleMatrix SSIq;
    public DoubleMatrix IqSS;
    public DoubleMatrix IqSI;
    public DoubleMatrix ISIq;
   
    public boolean SIR = false;
    public boolean estimateR0 = false;
    
    private double origin;
    private double timeStep;
    	
	@Override
	public void initAndValidate() throws Exception {
		
		states = maxKInput.get().getValue() + 1; 
		timeStep = timeStepInput.get().getValue();
		
		if(originInput.get() != null) origin  = originInput.get().getValue();
		
		/*
		 * Get the population model used as Input
		 */
		if (isSIRInput.get() != null) SIR = isSIRInput.get().getValue();
		
		/*
		 * If R0 is estimated than tau will be computed from R0
		 */
		if (estimateR0Input.get() != null) estimateR0 = estimateR0Input.get();

		// timeS is the timeSeries (refactor as timeSeries??)
		timeS = new TimeSeries();
	}
	
	public boolean update() throws MathException{
		
		timeS.setBack(states);	// empties all arrays in timeSeries
		
    	boolean reject = false;
//        if (!dirty) return reject;        
        
        final double endTime = 0.0; // time = 0 is the present
        
        // Pairwise model parameters
        final double nodes = nodesInput.get().getValue();
        tau = tauInput.get().getValue();
        nu = nuInput.get().getValue();
        R0 = R0Input.get().getValue();
        
        // Now defined through degreeDistInput
        //final double meanK = meanKInput.get().getValue();
        //final double varK = varKInput.get().getValue();
        
        // Get degree distribution pk
        pk = new DoubleMatrix(states); pk.put(0,0); // degree frequencies
        kSeries = new DoubleMatrix(states); kSeries.put(0,0); // degree series
        qSeries = new DoubleMatrix(states); qSeries.put(0,0); // only used for triple closures
        ParametricDistribution degreeDist = degreeDistInput.get();
        for (int k = 1; k < states; k++) {
        	final double kDouble = (double) k;
        	kSeries.put(k,kDouble);
        	pk.put(k, degreeDist.cumulativeProbability(kDouble) - degreeDist.cumulativeProbability(kSeries.get(k-1)));
        	qSeries.put(k, (kDouble - 1) / kDouble);
        }
        
        // Check for numerical stability
        //double pkMassThreshold = 0.95;
		//if (pk.sum() < pkMassThreshold) {
			//System.out.println("pk mass threshold exceded");
			//reject = true;
			//return reject;
		//}
        pk = pk.div(pk.sum()); // make sure pk sums to one
        //System.out.println(pk);
        
        // Compute moments of pk
        meanK = pk.dot(kSeries); // mean degree
        double pairs = meanK * nodes / 2; // number of pairs in net
        DoubleMatrix devs = kSeries.sub(meanK);
        DoubleMatrix sqrDevs = devs.mul(devs); // squared deviations
        varK = pk.dot(sqrDevs); // variance in pk
        
        // Get edge distribution ekl
        DoubleMatrix wpk = pk.mul(kSeries); // pk weighted by # of contacts
        DoubleMatrix ekl = wpk.mmul(wpk.transpose()); // outer product
        ekl.copy(ekl.div(meanK*meanK));
        
        // Compute transmission rate tau from R0
        if (estimateR0) {
        	DoubleMatrix Nkl = ekl.mul(pairs);
        	DoubleMatrix Nk = pk.mul(nodes);
        	tau = computeTauFromR(Nkl,Nk);
        }
        
        // Get initial node and pair state variables
        double expectedSI = meanK;
        double expectedSS = 2 * (pairs - expectedSI);
        DoubleMatrix pairsSS = ekl.mul(expectedSS);
        DoubleMatrix pairsSI = ekl.mul(expectedSI);
        DoubleMatrix pairsIS = pairsSI.transpose();
        DoubleMatrix pairsII = ekl.mul(0);
        
        // Set time series at t = 0
        timeS.setS(pk.mul(nodes).sub(pk));
        timeS.setI(pk);
        timeS.setSSPairs(pairsSS);
        timeS.setSIPairs(pairsSI);
        timeS.setISPairs(pairsIS);
        timeS.setIIPairs(pairsII);
        
        // Get integration times
        integrationTimes = new ArrayList<Double>();    	
		double time = origin;
		while (time > endTime){
			integrationTimes.add(time);
			time -= timeStep;
		}		
		integrationTimes.add(endTime);
	
		// Add to state variable trajectories
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
	 * Update all individual and pair-level state variables by Euler integration
	 */
	protected boolean deltaStep(int t, int dtTimeCount, double TimeStep){
		
		// Changes in pair variables
		DoubleMatrix dSS = new DoubleMatrix(states,states);
		DoubleMatrix dSI = new DoubleMatrix(states,states);
		DoubleMatrix dIS = new DoubleMatrix(states,states); // not needed
		DoubleMatrix dII = new DoubleMatrix(states,states);
		
		// Compute SSIq, IqSS, IqSI and ISIq
		approximateTriples();	
		
		// dSS				
		if (SIR) {
			dSS.copy(dSS.sub(SSIq.add(IqSS).mul(tau))); // transmission from SSI and ISS triples
		} else { // SIS
			dSS.copy(timeS.getSI().add(timeS.getIS()).mul(nu)); // add recoveries from SI and IS
			dSS.copy(dSS.sub(SSIq.add(IqSS).mul(tau))); // transmission from SSI and ISS triples
		}
		//DoubleMatrix dSSDelta = dSS.mul(tempDt);
		
		//dSI
		if (SIR) {
			dSI.copy(dSI.sub(timeS.getSI()).mul(nu));
			dSI.copy(dSI.add(SSIq.sub(IqSI).sub(timeS.getSI()).mul(tau)));
		} else { // SIS
			dSI.copy(timeS.getII().sub(timeS.getSI()).mul(nu));
			dSI.copy(dSI.add(SSIq.sub(IqSI).sub(timeS.getSI()).mul(tau)));
		}
		//DoubleMatrix dSIDelta = dSI.mul(tempDt);
		
		dII.copy(timeS.getII().mul(-2*nu));
		dII.copy(dII.add(ISIq.add(IqSI).add(timeS.getSI()).add(timeS.getIS()).mul(tau)));
		//DoubleMatrix dIIDelta = dII.mul(tempDt);
		
		// Changes in single variables
		DoubleMatrix dS = new DoubleMatrix(states);
		DoubleMatrix dI = new DoubleMatrix(states);
		DoubleMatrix dR = new DoubleMatrix(states);
		
		// Old way using columnSums (also works)
		//DoubleMatrix columnSumsSI = timeS.getSI().columnSums();
		//dS.copy(timeS.getI().mul(nu).sub(columnSumsSI.mul(tau)));
		//dI.copy(columnSumsSI.mul(tau).sub(timeS.getI().mul(nu)));
		
		// Doesn't matter if use rowSums or columnSums here
		DoubleMatrix rowSumsSI = timeS.getSI().rowSums();
		if (SIR) {
			dS.copy(dS.sub(rowSumsSI.mul(tau)));
		} else { // SIS
			dS.copy(timeS.getI().mul(nu).sub(rowSumsSI.mul(tau)));	
		}
		dI.copy(rowSumsSI.mul(tau).sub(timeS.getI().mul(nu)));
		
		// Transmission rates for coalescent model
		DoubleMatrix Transmission = timeS.getIS().mul(tau); // Need IS, not SI, so transmission is from k --> l
		
		// Have to add F, G and Y
		timeS.addF(Transmission);
		
		// No migration so no G matrix
		//DoubleMatrix Migration = migrationRate.mulColumnVector(timeS.getS());
		//timeS.addG(Migration);
		
		// Population sizes for coalescent model
		timeS.addY(timeS.getI());
		
		if (t < (dtTimeCount - 1)) {
		    
	    	// Update state variables
	    	final double dt = integrationTimes.get(t) - integrationTimes.get(t+1);
	    		
	    	//DoubleMatrix newS = timeS.getS().add(dS.mmul(tempDt));
	    	timeS.setS(timeS.getS().add(dS.mul(dt)));
	    	//DoubleMatrix newI = timeS.getI().add(dI.mmul(tempDt));
	    	timeS.setI(timeS.getI().add(dI.mul(dt)));
	    	
	    	//DoubleMatrix newSS = timeS.getSS().add(dSS.mul(dt)); // for debugging
	    	timeS.setSSPairs(timeS.getSS().add(dSS.mul(dt)));
	    	//DoubleMatrix newSI = timeS.getSI().add(dSI.mul(dt));  // for debugging
	    	timeS.setSIPairs(timeS.getSI().add(dSI.mul(dt)));
	    	//DoubleMatrix newIS = timeS.getSI().transpose();  // for debugging
	    	timeS.setISPairs(timeS.getSI().transpose());
	    	//DoubleMatrix newII = timeS.getII().add(dII.mul(dt));  // for debugging
	    	timeS.setIIPairs(timeS.getII().add(dII.mul(dt)));
	    	
	    	
	    	//timeS.setR(timeS.getR().add(dR.mmul(dt))); not track R
	    	//timeS.setN(timeS.getS().add(timeS.getI().add(timeS.getR()))); // constant	    			
	    	if (timeS.getPS().min() < 0|| timeS.getPI().min() < 0){
	    		return true;
	    	}
           	timeS.addS(timeS.getS());
           	timeS.addI(timeS.getI());
           	//timeS.addR(timeS.getR()); // not tracking R
           	
           	// Add add methods for pair variables
           	timeS.addSSPairs(timeS.getSS());
           	timeS.addSIPairs(timeS.getSI());
           	timeS.addISPairs(timeS.getIS());
           	timeS.addIIPairs(timeS.getII());
           	
           	/*
           	 * The following is only for debugging
           	 */
           	
           	//double pairCount = timeS.getSS().sum() + timeS.getSI().sum() + timeS.getIS().sum() + timeS.getII().sum(); 
           	//double indvCount = timeS.getS().sum() + timeS.getI().sum();
     
           	//DoubleMatrix currI = timeS.getI();
           	//DoubleMatrix currS = timeS.getS();
           	
           	//System.out.println(currI);
           	//System.out.println(currS);
           	//System.out.println(pairCount);
           	//System.out.println(indvCount);
           	
	    }
		
		return false;
	}
	
	
	/*
	 * Update triple approximations using moment closure (ignores clustering)
	 */
	protected void approximateTriples() {
		
		//These now all appear correct
		
		// With loops - right triple SSIq
		SSIq = new DoubleMatrix(states,states);
		DoubleMatrix cSS = timeS.getSS();
		DoubleMatrix rowSumsSI = timeS.getSI().rowSums(); // was rowSums
		for (int k = 1; k < states; k++) {
			for (int l = 1; l < states; l++) {
				final double v = qSeries.get(l) * cSS.get(k,l) * rowSumsSI.get(l) / timeS.getS().get(l); 
				SSIq.put(k,l,v);
			}
		}
		
		// With matrix-vector ops
		//DoubleMatrix SSI2 = cSS.mulRowVector(qSeries).mulRowVector(rowSumsSI).divRowVector(timeS.getS()); // works but end with NaNs if denominator is zero
		
		// Left triple IqSS
		IqSS = new DoubleMatrix(states,states);
		DoubleMatrix columnSumsIS = timeS.getIS().columnSums(); // was columnSums
		for (int k = 1; k < states; k++) {
			for (int l = 1; l < states; l++) {
				final double v = qSeries.get(k) * cSS.get(k,l) * columnSumsIS.get(k) / timeS.getS().get(k); 
				IqSS.put(k,l,v);
			}
		}
		
		// Left triple IqSI
		IqSI = new DoubleMatrix(states,states);
		DoubleMatrix cSI = timeS.getSI();
		for (int k = 1; k < states; k++) {
			for (int l = 1; l < states; l++) {
				final double v = qSeries.get(k) * cSI.get(k,l) * columnSumsIS.get(k) / timeS.getS().get(k); 
				IqSI.put(k,l,v);
			}
		}
		
		// Right triple ISIq
		ISIq = new DoubleMatrix(states,states);
		DoubleMatrix cIS = timeS.getIS();
		for (int k = 1; k < states; k++) {
			for (int l = 1; l < states; l++) {
				final double v = qSeries.get(l) * cIS.get(k,l) * rowSumsSI.get(l) / timeS.getS().get(l); 
				ISIq.put(k,l,v);
			}
		}
		
	}
	
	/*
	 * Solve for transmission rate tau using eigen decomposition of the next generation contact matrix
	 */
	protected double computeTauFromR(DoubleMatrix Nkl, DoubleMatrix Nk) {
		
		double tRate;
		DoubleMatrix M = new DoubleMatrix(states,states); // next generation contact matrix
		for (int k = 0; k < states; k++) {
			for (int l = 0; l < states; l++) {
				if (Nk.get(l) > 0) {
					final double val = Nkl.get(k,l) * qSeries.get(l) / Nk.get(l); 
					M.put(k,l,val);
				} else {
					M.put(k,l,0);
				}
			}
		}
		
		DoubleMatrix eigVals = Eigen.eigenvalues(M).real();
		final double lv = eigVals.max(); // leading eigenvalue
		tRate = R0 * nu / (lv - 1);
		
		return tRate;
		
	}
	
	
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
    	return timeS.getY(t).get(k);
    }
    
    public DoubleMatrix getY(int t) {
    	return timeS.getY(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in F
     * @param l column in F
    */
    public double getF(int t, int k, int l) {
    	return timeS.getF(t).get(k,l);
    }
    
    public DoubleMatrix getF(int t) {
    	return timeS.getF(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in G
     * @param l column in G
    */
    public double getG(int t, int k, int l) {
    	return timeS.getG(t).get(k,l);
    }
    
    public DoubleMatrix getG(int t) {
		return timeS.getG(t);
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

        out.print("meanK\t");
        out.print("varK\t");
        //out.print("tau\t");

    }

    public void log(int nSample, PrintStream out) {
    	
        out.format("%g\t", meanK);
    	out.format("%g\t", varK);

    }

    public void close(PrintStream out) {
    }
    
	
	
}

