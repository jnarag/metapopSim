import java.util.ArrayList;
import java.util.List;

public class extractGenealogy {

    int n_lineages = 100;



    public void sampleLineages(){

        //sample lineages from the simulation from which we will construct a transmission tree

    }

    public List<List<Integer>> getParentLineages(List<Integer> sampledLineages, List<Integer> parents){

        //get parents of particular set of lineages of interest

        List<List<Integer>> parentLineages = new ArrayList<List<Integer>>();

        for(int i=0; i<sampledLineages.size(); i++) {

            List<Integer> lineages = getSampleLineage(sampledLineages.get(i), parents);
            parentLineages.add(lineages);

        }

        return parentLineages;

    }

    public List<Integer> getSampleLineage(Integer indiv_lineage, List<Integer> parents){

        int i = 0;

        List<Integer> this_lineage = new ArrayList<Integer>();

        this_lineage.add(indiv_lineage);

        boolean condition = true;
        while(condition) {

            Integer lineage_parent = parents.get(this_lineage.get(i));
            if(lineage_parent < 0) {

                return this_lineage;

            }
            else{

                this_lineage.add(lineage_parent);
                i++;
            }
        }

        return this_lineage;

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

    public static void main(String [] args) {

        onePatch_test test = new onePatch_test();
        test.runSim();
        infectionHistory infectionHistory =  test.getInfectionHistory();

        System.out.println(infectionHistory.birth.get(infectionHistory.birth.size()-1));

        double moreRecentBirth = infectionHistory.birth.get(infectionHistory.birth.size()-1);

        List<List<Integer>> prevalenceHistory = infectionHistory.prevalence;

        int count = 1;
        for(int i=0; i < infectionHistory.getCompleteBirths().size(); i++) {

            if (infectionHistory.birth.get(i) > moreRecentBirth-50) {

                List<Integer> prevalence = prevalenceHistory.get(i);

                int lastPrev = -1;
                for(int p = 0; p < prevalence.size(); p++) {

                    if( prevalence.get(p) > 0) {

                        lastPrev = prevalence.get(p);
                    }

                }

                System.out.println(count+ " " +infectionHistory.genotype.get(i)
                                    + " " + (double)prevalenceHistory.get(i).lastIndexOf(lastPrev)/2.0
                                    + " " + infectionHistory.birth.get(i)
                                    + " " + infectionHistory.parent.get(i)
                                    + " " + prevalenceHistory.get(i));
                count++;
            }
        }


    }

}
