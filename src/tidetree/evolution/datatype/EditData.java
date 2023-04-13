package tidetree.evolution.datatype;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType.Base;

/**
 * @author Sophie Seidel
 */
@Description("Datatype for integer sequence representing the barcode states. The state 0 is unedited, " +
        "states 1, ..., N-1 represent particular indels, state N represents the silenced state.")
public class EditData extends Base {

    public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates",
            "The number of states a barcode can be in, i.e. the number of different indels + 2");


    public void initAndValidate(){
        stateCount  = nrOfStatesInput.get();
        mapCodeToStateSet = null;
        codeLength = 1;
        codeMap = null;
    }

    @Override
    public String getTypeDescription() {
        return "editData";
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
