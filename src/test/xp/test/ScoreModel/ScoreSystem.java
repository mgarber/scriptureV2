package xp.test.ScoreModel;

import nextgen.core.annotation.Annotation;

/**
 *  Created on 2013-3-5  
 */
public interface ScoreSystem {

	public Double getScore(Annotation a);
	public Double getScore(String chr,int start,int end);
	
}
