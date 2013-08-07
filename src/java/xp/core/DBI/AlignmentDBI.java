package xp.core.DBI;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 *  Created on 2013-3-4  
 */
public interface AlignmentDBI<T> {
	
   Iterator<T> query(String chr,int start ,int stop);
   Iterator<T> iterate();
   int getCount(String chr,int start,int stop);
   double getGlobalCount();
   double getLocalCount(String chr,int start,int stop);
   List<String>  getChrs();
   long  getChrLength(String chr);
    
}


