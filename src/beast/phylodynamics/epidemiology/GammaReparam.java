package beast.phylodynamics.epidemiology;


import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;



@Description("Reparameterization of the Gamma distribution. Alpha and beta are now compute based on input mean and variance."
		+ " for x>0  g(x;alpha,beta) = \\frac{beta^{alpha}}{Gamma(alpha)} x^{alpha-1} e^{-beta {x}}" +
        "If the input x is a multidimensional parameter, each of the dimensions is considered as a " +
        "separate independent component.")
public class GammaReparam extends ParametricDistribution {
    
	public Input<RealParameter> meanInput = new Input<RealParameter>("mean", "shape parameter, defaults to 2");
    public Input<RealParameter> varInput = new Input<RealParameter>("variance", "scale parameter, defaults to 2");

    static org.apache.commons.math.distribution.GammaDistribution m_dist = new GammaDistributionImpl(1, 1);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    void refresh() {
    	
        double fMean;
        double fVar;
        if (meanInput.get() == null) {
            fMean = 2;
        } else {
            fMean = meanInput.get().getValue();
        }
        if (varInput.get() == null) {
            fVar = 2;
        } else {
            fVar = varInput.get().getValue();
        }
        
        final double beta = fVar / fMean;
        final double alpha = fMean / beta; 
        
        m_dist.setAlpha(alpha);
        m_dist.setBeta(beta);
    }

    @Override
    public ContinuousDistribution getDistribution() {
        refresh();
        return m_dist;
    }

    @Override
    public double getMean() {
    	return offsetInput.get() + m_dist.getAlpha() / m_dist.getBeta();
    }
} // class Gamma
