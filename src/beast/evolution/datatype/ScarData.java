package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType.Base;

@Description("Datatype for integer sequence representing sequence states.")
public class ScarData extends Base {

    public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates",
            "The number of states a sequence can be in.");


    public void initAndValidate(){
        stateCount  = nrOfStatesInput.get();
        mapCodeToStateSet = null;
        codeLength = 1;
        codeMap = null;
        System.out.print("This is the number for scarDATA" + codeLength + "\n");
    }

    @Override
    public String getTypeDescription() {
        return "scarData";
    }

    @Override
    public boolean isAmbiguousCode(int code) {
        return code < 0;
    }

    @Override
    public String getCharacter(int code) {
        if (code < 0) {
            return "?";
        }
        return code + "";
    }

    @Override
    public int[] getStatesForCode(int code) {

        return new int[]{code};
    }
}
