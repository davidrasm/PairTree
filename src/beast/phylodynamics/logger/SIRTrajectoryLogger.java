package beast.phylodynamics.logger;

import java.io.PrintStream;

import org.apache.commons.math.MathException;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.phylodynamics.model.PairwiseEpiModel;


/**
 * @author Nicola Mueller / David Rasmussen
 */
@Description("Logs the SIR variable trajectories")
public class SIRTrajectoryLogger extends BEASTObject implements Loggable {
    
	public Input<PairwiseEpiModel> epiModelInput = new Input<PairwiseEpiModel>(
    		"epiModel",
            "the rates of migration between and within demes, the dimension needs to be n x n.",
            Validate.REQUIRED);
    public Input<RealParameter> printEveryInput = new Input<RealParameter>(
    		"printEvery",
    		"how strong the seasonal patterns are, higher means less seasonal"); 
    
    protected int states;
    private int printEvery;
    
    @Override
    public void initAndValidate() throws Exception {
    	double tmp;
    	if (printEveryInput.get()!=null) tmp=printEveryInput.get().getValue();
    	else tmp=1.0;
    	printEvery = (int) tmp;
    }
	

    @Override
    public void init(PrintStream out) throws Exception {
    	epiModelInput.get().update();
    	for (int i = 0; i < epiModelInput.get().getTrajLength(); i=i+printEvery){
			if ( i < epiModelInput.get().getTrajLength()-1){
				out.print("t_" + i + "\t");
			}else{
				out.print("t_" + i);
			}
    	}
    	out.print("\n");
    	out.print("0" + "\t");
    	out.print("0" + "\t");
    	for (int i = printEvery; i < epiModelInput.get().getTrajLength(); i=i+printEvery){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.print(epiModelInput.get().getTime(i) + "\t");
				}else{
					out.print(epiModelInput.get().getTime(i));
				}
	    	}
    }
    
    @Override
    public void log(int nSample, PrintStream out) {
    	try {
			epiModelInput.get().update();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

        //Susceptible Fraction;
//    	for (int d = 0; d < epiModelInput.get().states; d++){
//    		if (d==0){
//    		}
//    		else{
//    			out.print(nSample + "\t");
//    		}
//    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i=i+printEvery){
//				if ( i < epiModelInput.get().getTrajLength()-1){
//					out.format("%.2f\t",epiModelInput.get().getNumSDeme(i, d));
//				}else{
//					out.format("%.2f",epiModelInput.get().getNumSDeme(i, d));
//				}
//	    	}
//	    	out.print("\n");
//    	}
    	
    	//Infected Fraction;
    	//for (int d = 0; d < epiModelInput.get().states; d++){
    		//if (d==0){
    		//}
    		//else{
    		//	out.print(nSample + "\t");
    		//}
    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i=i+printEvery){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.format("%.2f\t",epiModelInput.get().getNumI(i).sum());
					//out.format("%.2f\t",epiModelInput.get().getNumIDeme(i, d));
				}else{
					out.format("%.2f",epiModelInput.get().getNumI(i).sum());
					//out.format("%.2f",epiModelInput.get().getNumIDeme(i, d));
				}				
	    	}
	    	out.print("\n");
    	//}
    	
//    	//F matrix;
//    	for (int d = 0; d < epiModelInput.get().states; d++){
//   			out.print(nSample + "\t");
//    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i=i+printEvery){
//				if ( i < epiModelInput.get().getTrajLength()-1){
//					out.format("%.2f\t",epiModelInput.get().getF(i,d, d));
//				}else{
//					out.format("%.2f",epiModelInput.get().getF(i,d, d));
//				}		
//	    	}
//    		if(d < (epiModelInput.get().states-1)){
//    			out.print("\n");
//    		}
//    	}
    }
    
    @Override
    public void close(PrintStream out) {
    }

}
