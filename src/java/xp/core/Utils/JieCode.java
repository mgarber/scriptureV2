package xp.test.Utils;
/**
 *  Created on 2013-3-7  
 */
public class JieCode {
	private int tid; //chromosome tid2chr
	private int pos; //0-index
	private boolean isStart;  //isStart or End
	private int classIndex;    // belong to which class;
	private int USB;   // for future using
	private static int DEFAULT_USB=0;
	public JieCode(int tid, int pos, boolean isStart, int classIndex) {
		super();
		this.tid = tid;
		this.pos = pos;
		this.isStart = isStart;
		this.classIndex = classIndex;
		this.USB=DEFAULT_USB;
	}
	public JieCode(int tid, int pos, boolean isStart, int classIndex, int USB) {
		super();
		this.tid = tid;
		this.pos = pos;
		this.isStart = isStart;
		this.classIndex = classIndex;
		this.USB = USB;
	}
	public int getTid() {
		return tid;
	}
	public void setTid(int tid) {
		this.tid = tid;
	}
	public int getPos() {
		return pos;
	}
	public void setPos(int pos) {
		this.pos = pos;
	}
	public boolean isStart() {
		return isStart;
	}
	public void setStart(boolean isStart) {
		this.isStart = isStart;
	}
	public int getClassIndex() {
		return classIndex;
	}
	public void setClassIndex(int classIndex) {
		this.classIndex = classIndex;
	}
	public int getUSB() {
		return USB;
	}
	public void setUSB(int USB) {
		this.USB = USB;
	}
	

}
