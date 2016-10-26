package beast.phylodynamics.epidemiology;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;
import beast.math.distributions.ParametricDistribution;
import beast.phylodynamics.model.PairwiseEpiModel;


/**
 * @author David Rasmussen
 */

@Description("Calculates the probability of a beast.tree using a pairwise network-based epidemic/coalescent model.")
public class StructCoalNetDensity extends StructuredTreeDistribution implements Loggable {

	public Input<PairwiseEpiModel> epiModelInput = new Input<PairwiseEpiModel>(
			"epiModel",
			"A deterministic pairwise epidemiological model",
			Input.Validate.REQUIRED);
	
    public Input<TraitSet> typeTraitInput = new Input<TraitSet>(
            "typeTrait", "Type trait set.");
    
    public Input<BooleanParameter> tipStatesProvided = new Input<BooleanParameter>(
    		"tipStates",
    		"define whether tip states are provided");
    
    public Input<BooleanParameter> approxLambdaSumInput = new Input<BooleanParameter>(
    		"approxLambdaSum",
    		"define whether to compute lambda sum by binning lineages by state");
    
    public Input<BooleanParameter> ghostLineCorrectionInput = new Input<BooleanParameter>(
    		"ghostLineCorrection",
    		"define whether to use ghost lineage correction when updating lineage state probabilities.");
    
    public Input<RealParameter> minIkIlInput = new Input<RealParameter>(
    		"minIkIl",
    		"Set minimum threshold for IkIl denominator in coal rate to smooth likelihood.");
    
    public Input<Boolean> fixOriginAtRootInput = new Input<Boolean>(
            "fixOriginAtRoot",
            "Is the newick tree labelled? Default=false.", false);
	
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
	
	public DoubleMatrix tipPk;
    
    public int states;
    
    private boolean needsUpdate = true;
    private boolean traitInput = false;
    private boolean tipStates = false;
    private boolean ghostLineCorrection = false;
    private boolean approxLambdaSum = false;
    private double minIkIl;
    private boolean fixOriginAtRoot = false;
    
    // Set up for lineage state probabilities
    ArrayList<Integer> activeLineages;
    ArrayList<DoubleMatrix> lineStateProbs;
    
    @Override
    public void initAndValidate() throws Exception {
    	
    	treeIntervalsInput.get().calculateIntervals();
        if (treeIntervalsInput.get() == null)
            throw new Exception("Expected treeIntervals to be specified");
        
        stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];               
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;

        states = epiModelInput.get().states;
        
        if (typeTraitInput.get() != null) traitInput = true;
        if (tipStatesProvided.get() != null) tipStates = tipStatesProvided.get().getValue();
        if (approxLambdaSumInput.get() != null) approxLambdaSum = approxLambdaSumInput.get().getValue();
        if (ghostLineCorrectionInput.get() != null) ghostLineCorrection = ghostLineCorrectionInput.get().getValue();
        if (minIkIlInput.get() != null) minIkIl = minIkIlInput.get().getValue();
        if (fixOriginAtRootInput.get() != null) fixOriginAtRoot = fixOriginAtRootInput.get();
           
        calculateLogP();
    }

    public double calculateLogP() throws Exception {

    	boolean reject = epiModelInput.get().update();
    	treeIntervalsInput.get().requiresRecalculation();
    	
    	if (fixOriginAtRoot) {
    		
    		epiModelInput.get().setOrigin(treeIntervalsInput.get().getTotalDuration());
    		
    	} else {
    	
	    	// if origin is being sampled it must be constrained to be older than the root of the tree.
	    	if (epiModelInput.get().getOrigin() < treeIntervalsInput.get().getTotalDuration()) {
		    	System.out.println("origin is too low");
		        logP = Double.NEGATIVE_INFINITY;
		        System.out.println("Coal logP is NEG INF");
		        return logP;
	        }
    	
    	}
    	
        // if maximal count exceeded or number too small
        if (reject) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }
        
        // Set up for lineage state probabilities
        activeLineages = new ArrayList<Integer>();
        lineStateProbs = new ArrayList<DoubleMatrix>();
        
        // Compute likelihood at each integration time and tree event starting at final sampling time and moving backwards
        ArrayList<Double> integrationTimes = epiModelInput.get().getIntegrationTimes();							//in reverse time
        int currTreeInterval = 0; 																				// what tree interval are we in?
        double currTime = integrationTimes.get(0); 																// current time (since present)
        double nextIntegrationTime; 																			// time of next integration step
        double nextIntervalTime = treeIntervalsInput.get().getInterval(0); 										// time next tree interval starts
        final int intervalCount = treeIntervalsInput.get().getIntervalCount();
        
        logP = 0;
        int t = 0;
        
        do {        	
        	nextIntegrationTime = integrationTimes.get(t+1);													// Update Integration time
        	while (nextIntervalTime <= nextIntegrationTime) {													// while there are still events to take place before the next integration step
        		/*
        		 * compute contribution of last interval to likelihood
        		 */
	        	final double duration = nextIntervalTime - currTime;
	        	if (duration > 0) {
	        		double add = -computeLambdaSum(t) * duration; // this is still done the slow way
	        		if (add > 0)
	        		{
	        			System.out.println("lambda sum is larger than one:\t" + add + "\t" + duration);
	        			epiModelInput.get().getF(t).print();
	        			epiModelInput.get().getY(t).print();
	        		}	
	        		logP += add;
	        		
	        	}                
	               
	        	/*
	        	 *  compute contribution of event to likelihood
	        	 */
	        	if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {	//CHECK IF INDEX IS CORRECT HERE
	        		
	        		final double yTotal = computeYTotal(t);
	        		if (yTotal < 1.0) {
	        			logP = Double.NEGATIVE_INFINITY;
	        			return logP;
	        		}
	        		
	        		final double pairCoalRate = computeLambdaPair(t, currTreeInterval);
	        		if (pairCoalRate >= 0) {
	        			logP += Math.log(pairCoalRate);
	        		}
	        		updateLineProbsCoalEvent(t, currTreeInterval);												// Set parent lineage state probs and remove children
		       	}
	       		
	        	/*
	        	 * add new lineage
	        	 */
	       		if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.SAMPLE) {
	       			addLineages(t,currTreeInterval); // added t input so can initialize line states based on pop vars
	       		}
	       		
	       		currTime += duration;
	       		currTreeInterval++; 
	       		
	       		if (currTreeInterval < intervalCount) {
	       			nextIntervalTime = currTime + treeIntervalsInput.get().getInterval(currTreeInterval);
        		} else {
	       			nextIntervalTime = Double.POSITIVE_INFINITY;
	       			currTreeInterval--; //stay in last interval
	       		}

        	}   
        	
	        final double duration = nextIntegrationTime - currTime;
	    	if (duration > 0) {    			
	   			// compute contribution of last time step to likelihood
  				logP += -computeLambdaSum(t) * duration; // if A is greater than Y???
  				boolean negInf = updateLineProbs(t, duration); 
  				if (negInf)
  					logP = Double.NEGATIVE_INFINITY;
    		}		    	
    		currTime += duration;
    		t++;
    		
        }while((integrationTimes.get(t-1) <= treeIntervalsInput.get().getTotalDuration()) && (t<integrationTimes.size()-1)); // added second condition for case where root if final integration time
	        
        //Compute likelihood of remaining tree intervals (coal events occuring before origin)
        if (Double.isInfinite(logP))logP = Double.NEGATIVE_INFINITY;
        return logP;
   	
    }
    
	/*
	 * Compute the pairwise coalescent rate for lineage pair
	 */
    private double computeLambdaPair(int t, int currTreeInterval) throws Exception {
		
		List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			throw new Exception("Unsupported coalescent at non-binary node");
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			return 0.0;
		}
		double lambda = 0.0;
		
        for (int k = 0; k < states; k++) { 
        	double Yk = epiModelInput.get().getY(t,k);
            if ( Yk > 0.0) {
            	for (int l = 0; l < states; l++){
                	double Yl = epiModelInput.get().getY(t,l);
                	if (Yl > 0.0){
                		
						// Modified for PairNetCoal models
                		double pairCoalRate;
						if (k == l) {
						
							double iChoose2;
							if (Yk > 1.0) {
								iChoose2 = Yk * (Yk - 1) / 2;
							} else {
								iChoose2 = Yk * Yk / 2; // don't need this if taking max
							}
							
							iChoose2 = Math.max(minIkIl,iChoose2);
							final double popCoalRate = epiModelInput.get().getF(t,k,l) / iChoose2;
							
							pairCoalRate = popCoalRate * (lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l));
									//+ lineStateProbs.get(daughterIndex1).get(l) * lineStateProbs.get(daughterIndex2).get(k)); // modified this to check pair coal rate for k = l
							
						} else {
							
							final double IkIl = Math.max(minIkIl,Yk*Yl);
							final double popCoalRate = epiModelInput.get().getF(t,k,l) / IkIl;
							pairCoalRate = popCoalRate * (lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l)
									+ lineStateProbs.get(daughterIndex1).get(l) * lineStateProbs.get(daughterIndex2).get(k));
							
						}
						if (!Double.isNaN(pairCoalRate)) lambda += pairCoalRate;
						
                	}
            	}
            }
        }
		
    	return lambda;
    }
    
	/*
	 * Compute the coalescent rate for all lineages
	 */
    private double computeLambdaSum(int t) {
    	
		int lineCount = activeLineages.size();
		double lambdaSum = 0.0; 										// holds the sum of the pairwise coalescent rates over all lineage pairs
		
		if (lineCount < 2) return lambdaSum;											// Zero probability of two lineages coalescing	


		/*
		 * Sum over line pairs (scales better with numbers of demes and lineages)
		 */	
		DoubleMatrix popCoalRates;
		
		/*
		 * If approxLambdaSum, approximate lambda sum over lineage pairs by binning lineages by state
		 */
    	if(approxLambdaSum){
    		
    		DoubleMatrix A = getLineStateSum();
			for (int k = 0; k < states; k++) {
				for (int l = 0; l < states; l++) {
					final double Yk = epiModelInput.get().getY(t,k);
					final double Yl = epiModelInput.get().getY(t,l);
					if (Yk > 0 & Yl > 0) {
						
						// Modified for PairNetCoal models
						if (k == l) {
							
							double pairsAkAk; // = A.get(k) * A.get(l) / 2;
							if (A.get(k) > 1) {
								pairsAkAk = A.get(k) * (A.get(k) - 1) / 2;
							} else {
								pairsAkAk = A.get(k) * A.get(k) / 2;
							}
							
							double iChoose2;
							if (Yk > 1.0) {
								iChoose2 = Yk * (Yk - 1) / 2;
							} else {
								iChoose2 = Yk * Yk / 2; // don't need this if taking max
							}
							
							iChoose2 = Math.max(minIkIl,iChoose2);
							final double lambda = pairsAkAk * epiModelInput.get().getF(t,k,l) / iChoose2;
							if (!Double.isNaN(lambda)) lambdaSum += lambda;
							
						} else {
							
							final double pairsAkAl = A.get(k) * A.get(l); 
							final double IkIl = Math.max(minIkIl,Yk*Yl);
							final double lambda = pairsAkAl * epiModelInput.get().getF(t,k,l) / IkIl;
							if (!Double.isNaN(lambda)) lambdaSum += lambda;
							
						}
						
					}
				}	
			}
			
		/*
		 * else sum over all lineage pairs
		 */
    	}else{    		
			popCoalRates = new DoubleMatrix(states, states);
			for (int k = 0; k < states; k++) {
				for (int l = 0; l < states; l++) {
					final double Yk = epiModelInput.get().getY(t,k);
					final double Yl = epiModelInput.get().getY(t,l);
					if (Yk > 0 & Yl > 0) {
						
						// Modified for PairNetCoal models
						if (k == l) {
						
							double iChoose2;
							if (Yk > 1.0) {
								iChoose2 = Yk * (Yk - 1) / 2;
							} else {
								iChoose2 = Yk * Yk / 2;
							}
							iChoose2 = Math.max(minIkIl,iChoose2);
							popCoalRates.put(k,l, (epiModelInput.get().getF(t,k,l) / iChoose2));
							
						} else {
							
							final double IkIl = Math.max(minIkIl,Yk*Yl);
							popCoalRates.put(k,l, (epiModelInput.get().getF(t,k,l) / IkIl));
							
						}
						
						
					} else {
						popCoalRates.put(k,l, 0);
					}
				}	
			}
			
			for (int linI = 0; linI < lineCount; linI++) {					
				for (int linJ = linI+1; linJ < lineCount; linJ++) {

					
					//Sum over line states
					for (int k = 0; k < states; k++) {
						for (int l = 0; l < states; l++) {
							
							if (k==l){	
								final double pairCoalRate = popCoalRates.get(k,l) * (lineStateProbs.get(linI).get(k) * lineStateProbs.get(linJ).get(l));
								if (!Double.isNaN(pairCoalRate)) lambdaSum += pairCoalRate;
							} else {
								final double pairCoalRate = popCoalRates.get(k,l) * (lineStateProbs.get(linI).get(k) * lineStateProbs.get(linJ).get(l)
										+ lineStateProbs.get(linI).get(l) * lineStateProbs.get(linJ).get(k));
								if (!Double.isNaN(pairCoalRate)) lambdaSum += pairCoalRate;	
							}

						}
					}
				}	
			}
    	}
    	return lambdaSum;    	
    }
    
	/*
	 * Update lineage state probabilities
	 */
    private boolean updateLineProbs(int t, double dt) {    	

		DoubleMatrix A = getLineStateSum();
		DoubleMatrix Y = epiModelInput.get().getY(t);
		//DoubleMatrix mul = Y.sub(A).div(Y);  // not currently used now
		DoubleMatrix M = new DoubleMatrix();
			
		/*
		 * Set up the Q matrix to hold lineage transition rates
		 */
    	DoubleMatrix Q = DoubleMatrix.zeros(states,states);
   	
		for (int k = 0; k < states; k++) 
		{
			double rowSum = 0.0;
			for (int l = 0; l < states; l++) 
			{
				if (k != l) 
				{																        				// off-diagonal
					if (Y.get(k) > 0) 
					{
						
						double rateOutByBirth;
						if (ghostLineCorrection) {
							
							double pairProb;
							double corrIk;
							if (Y.get(l) > 0) {
								pairProb = 1 / (Y.get(k) * Y.get(l));
								corrIk = (Y.get(l) - A.get(l)) / Y.get(l);
								corrIk = Math.max(0,corrIk);
							} else {
								pairProb = 0;
								corrIk = 0;
							}
							rateOutByBirth = epiModelInput.get().getF(t,l,k) * pairProb * corrIk;       // Birth in other deme (Backwards in time)
						
						} else { 
							
							rateOutByBirth = epiModelInput.get().getF(t,l,k) / Y.get(k);       // Birth in other deme (Backwards in time)
							
						}
						//final double rateOutByMigration = epiModelInput.get().getG(t,l,k) / (Y.get(k)); 								// State change of lineages (Backwards in time)
						final double totalRate = rateOutByBirth; // + rateOutByMigration;
						Q.put(k, l, totalRate);
						rowSum += totalRate;
					} 
					else {																// diagonal
						Q.put(k, l, 0.0);
					}
				}
			}
			Q.put(k,k,-rowSum);
		}
		M = Q.mul(dt);    	

		/*
		 * Update lineage state probabilities for all active lineages. 
		 */
		DoubleMatrix newP = MatrixFunctions.expm(M);
		for (int lin = 0; lin < activeLineages.size(); lin++) {
			
			DoubleMatrix oldProbs = lineStateProbs.get(lin);
			DoubleMatrix probs = vecMat(newP, oldProbs);
			
			/*
			 * Vector-matrix mult gives similar probs
			 * But is actually slower when compared against matrix exp
			 * It may also be slightly less accurate but didn't test thoroughly
			 */
			//DoubleMatrix dProbs = oldProbs.transpose().mmul(Q).mul(dt);
			//DoubleMatrix probs = oldProbs.add(dProbs);
			//DoubleMatrix probDiffs = probs.sub(myProbs);
			//System.out.println(probDiffs);
			
			/*
			 * divide the state probabilities by the sum of state probabilities
			 * in order to avoid numerical errors
			 */
			probs = probs.div(probs.sum());
			lineStateProbs.set(lin, probs);
		}
		
		return false;
		
    }
    
	/*
	 * Update lineage state probabilities at a coalescent event for parent
	 */
    private void updateLineProbsCoalEvent(int t, int currTreeInterval) throws Exception {
    	
    	List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			throw new Exception("Unsupported coalescent at non-binary node");
		}
		int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		
		// The following ensures coalescent events occur in correct parent-offspring order if there are simultaneuous branching events in tree
		if (daughterIndex1 == -1 || daughterIndex2 == -1){
	    	treeIntervalsInput.get().swap();
	    	coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
	    	daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
	    	daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
	    	if (daughterIndex1 == -1 || daughterIndex2 == -1)
				throw new Exception("Active lineages does not contain coalescing lineages");
		}
		List<Node> parentLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if (parentLines.size() > 1) throw new Exception("Unsupported coalescent at non-binary node");
		
		//Add parent to activeLineage and initialize parent's state prob vector
		Node parentNode = parentLines.get(0);
		activeLineages.add(parentNode.getNr());
		DoubleMatrix pVec = new DoubleMatrix(states);

		//Compute parent lineage state probabilities in pVec
		DoubleMatrix coalRates = DoubleMatrix.zeros(states, states);
		double lambda;
		for (int k = 0; k < states; k++) {
            for (int l = 0; l < states; l++) {            	
            	final double Yk = epiModelInput.get().getY(t,k);
            	final double Yl = epiModelInput.get().getY(t,l);
            	if (Yk > 0 && Yl > 0) {
            		
					if (k == l) {
					
						double iChoose2;
						if (Yk > 1.0) {
							iChoose2 = Yk * (Yk - 1) / 2;
						} else {
							iChoose2 = Yk * Yk / 2;
						}
						// iChoose2 = Math.max(minIkIl,iChoose2);
						final double popCoalRate = epiModelInput.get().getF(t,k,l) / iChoose2;
						lambda = popCoalRate * (lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l));
						
					} else {
						
						//final double IkIl = Math.max(minIkIl,Yk*Yl);
						final double IkIl = Yk * Yl;
						final double popCoalRate = epiModelInput.get().getF(t,k,l) / IkIl;
						lambda = popCoalRate * (lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l)
								+ lineStateProbs.get(daughterIndex1).get(l) * lineStateProbs.get(daughterIndex2).get(k));
						
					}

            	} else { 
            		lambda = 0.0;
            	}
            	coalRates.put(k, l, lambda);
            }
        }	
		
		/*
		 * pVec := new state probabilities of the lineage dt after coalescing		
		 */
		pVec = coalRates.rowSums().div(coalRates.sum());
		stateProbabilities[parentNode.getNr() - nrSamples] = (pVec);
        lineStateProbs.add(pVec);
        
		//Remove daughter lineages
		activeLineages.remove(daughterIndex1);
		lineStateProbs.remove(daughterIndex1);
		
		daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr()); //index may have changed after removal of first daughter
		activeLineages.remove(daughterIndex2);
		lineStateProbs.remove(daughterIndex2);
    }
    
    private void addLineages(int t, int currTreeInterval) {
    	
		List<Node> incomingLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if(traitInput){
			/*
			 * If there is a typeTrait given as Input the model will take this
			 * trait as states for the taxons
			 */		
			for (Node l : incomingLines) {				
				activeLineages.add(l.getNr());
				//int sampleState = (int) typeTraitInput.get().getValue(l.getID()); // how do we fix this?
				int sampleState = 0;
				DoubleMatrix sVec = DoubleMatrix.zeros(states);
				sVec.put(sampleState, 1.0);			
				lineStateProbs.add(sVec);
			}			
		}else{		
			/*
			 * If there is no trait given as Input, the model will simply assume that
			 * the last value of the taxon name, the last value after a _, is an integer
			 * that gives the type of that taxon
			 */
			for (Node l : incomingLines) {
				activeLineages.add(l.getNr());
				
				DoubleMatrix sVec;
				if (tipStates) {
					String sampleID = l.getID();
					String[] splits = sampleID.split("-");
					int sampleState = Integer.parseInt(splits[0]); //samples states (or priors) should eventually be specified in the XML
					sVec = DoubleMatrix.zeros(states);
					sVec.put(sampleState, 1.0);
				} else {
					// If based on degree dist of infected pop Ik
					DoubleMatrix Yt = epiModelInput.get().getY(t);
					sVec = Yt.div(Yt.sum());
				}
				
				lineStateProbs.add(sVec);
			}	
		}
    }
    
    private DoubleMatrix getLineStateSum() {    	
    	DoubleMatrix sumAk = DoubleMatrix.zeros(states);
   		for (int lin = 0; lin < activeLineages.size(); lin++)
   			sumAk = sumAk.add(lineStateProbs.get(lin));
    	return sumAk;    	
    }
    
    private double computeYTotal(int t) {
    	return epiModelInput.get().getY(t).sum();
    }
    
    public DoubleMatrix getStateProb(int nr){
    	return stateProbabilities[nr - nrSamples];
    }
      
    
    public DoubleMatrix[] getStateProbabilities(){
    	return stateProbabilities;
    }
    
    public String getType(){
    	if (typeTraitInput.get()!=null) return typeTraitInput.get().getTraitName();
    	else return "type";
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	return true;
    }
    
    /*
     * Quickly added vector times matrix 
     */
    private DoubleMatrix vecMat(DoubleMatrix mat, DoubleMatrix vec){
    	DoubleMatrix prod = DoubleMatrix.zeros(states);
    	for (int i = 0; i < states; i++)
    		for (int j = 0; j < states; j++)
    			prod.put(i, prod.get(i) + vec.get(j)*mat.get(j,i));
    	return prod;
    }
    
}