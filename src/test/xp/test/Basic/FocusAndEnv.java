package xp.test.Basic;

import java.util.Arrays;


/**
 *  Created on 2013-3-9  
 */
public class FocusAndEnv {
   private BedGraphMultiScore focus;
   private BedGraphMultiScore envs;
public FocusAndEnv(BedGraphMultiScore focus, BedGraphMultiScore envs) {
	super();
	this.focus = focus;
	this.envs = envs;
}
@Override
public String toString() {
	return "FocusAndEnv [focus=" + focus + ", envs=" + envs + "]";
}
public BedGraphMultiScore getFocus() {
	return focus;
}
public void setFocus(BedGraphMultiScore focus) {
	this.focus = focus;
}
public BedGraphMultiScore getEnvs() {
	return envs;
}
public void setEnvs(BedGraphMultiScore envs) {
	this.envs = envs;
}
   
}