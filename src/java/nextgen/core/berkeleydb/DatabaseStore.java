package nextgen.core.berkeleydb;

import org.apache.log4j.Logger;

import com.sleepycat.je.Environment;
import com.sleepycat.persist.EntityStore;
import com.sleepycat.persist.StoreConfig;
import com.sleepycat.persist.evolve.Mutations;

/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_first.html
 */
/**
 * A database store
 * @author prussell
 *
 */
public class DatabaseStore {
	
	private EntityStore store;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(DatabaseStore.class.getName());
	
	public DatabaseStore() {}
	
	/**
	 * Set up the store
	 * @param env Database environment
	 * @param storeName Entity store name
	 * @param readOnly Whether the store should be read only
	 */
	public void setup(Environment env, String storeName, boolean readOnly) {
		setup(env, storeName, null, readOnly);
	}
	
	/**
	 * Set up the store
	 * @param env Database environment
	 * @param storeName Entity store name
	 * @param Mutations Collection of mutations for configuring class evolution
	 * @param readOnly Whether the store should be read only
	 */
	public void setup(Environment env, String storeName, Mutations mutations, boolean readOnly) {
		StoreConfig storeConfig = new StoreConfig();
		storeConfig.setReadOnly(readOnly);
		storeConfig.setAllowCreate(!readOnly);
		if(mutations != null) {
			storeConfig.setMutations(mutations);
		}
		store = new EntityStore(env, storeName, storeConfig);
	}
	
	/**
	 * Get the underlying EntityStore object
	 * @return
	 */
	public EntityStore getStore() {
		return store;
	}
	
	/**
	 * Close the store
	 */
	public void close() {
		if(store != null) {
			store.close();
		}
	}
	
}
