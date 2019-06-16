import java.util.List;

public class extractGenealogy {

    int n_lineages = 100;



    public void sampleLineages(){

    }

    public List<List<Integer>> getParentLineages(){

        return null;
    }

    public List<List<Integer>> getSampleLineages(){
        return null;

    }

    public coalescentEvent findMostRecentCoalescence() {

        return null;
    }

    class coalescentEvent{

        double timeToCoalescence;
        int coalParent;
        int[] coalDaughters;

        coalescentEvent(double timeCoalescence, Integer coalParent, int[] coalDaughters) {

            this.timeToCoalescence = timeCoalescence;
            this.coalParent = coalParent;
            this.coalDaughters = coalDaughters;
        }

        private double getTimeToCoalescence(){

            return this.timeToCoalescence;
        }
        private Integer getCoalParent() {

            return this.coalParent;
        }
        private int[] getCoalDaughters() {

            return this.coalDaughters;
        }
    }

}
