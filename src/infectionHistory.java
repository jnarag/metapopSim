import java.util.ArrayList;
import java.util.List;

/**
 * Created by jayna on 18/05/2018.
 */
public class infectionHistory {

    List<Integer> genotype = null;
    List<Integer> mutationsFromParent = null;
    List<Double> birth = null;
    List<Double> death = null;
    List<Integer> parent = null;
    List<List<Integer>> prevalence = null;
    List<Integer> patch = null;
    List<Double> fitness = null;

    public infectionHistory () {

        genotype = new ArrayList<>();
        mutationsFromParent = new ArrayList<>();
        birth = new ArrayList<>();
        death = new ArrayList<>();
        parent = new ArrayList<>();
        prevalence = new ArrayList<>();
        patch = new ArrayList<>();

    }

    public void logGenotype(Integer genotype) {

        this.genotype.add(genotype);

    }


    public void logMutations(Integer mutations) {

        this.mutationsFromParent.add(mutations);

    }

    public void logBirth(Double birth) {

        this.birth.add(birth);
    }

    public void logDeath(Double death) {

        this.death.add(death);
    }

    public void logParent(Integer parent) {

        this.parent.add(parent);
    }

    public void logPatch(Integer patch) {

        this.patch.add(patch);
    }




    public void logPrevalence(Integer prevalence, int index) {

        this.prevalence.get(index).add(prevalence);
    }
    public void logFitness(Double fitness) {

        this.fitness.add(fitness);

    }

    public int getId(int index) {

        return this.genotype.get(index);
    }

    public Double getBirth(int index) {

        return this.birth.get(index);
    }


    public Double getBirth(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.birth.get(index);
    }

    public Integer getMutationsFromParent(Integer genotype) {


        int index = this.genotype.indexOf(genotype);
        return this.mutationsFromParent.get(index);
    }

    public Integer getMutationsFromParent(Integer genotype, Integer parent) {

        if(genotype.equals(parent)) {
            return new Integer(0);
        }
        else {
            Integer mutations = getMutationsFromParent(genotype);
            Integer parent_i = getParent(genotype);


            while (!parent_i.equals(parent)) {

                mutations += getMutationsFromParent(parent_i);
                parent_i = new Integer(getParent(parent_i));

            }

            return mutations;
        }
    }

    public Double getDeath(int index) {

        return this.death.get(index);
    }

    public Integer getPatch(int index) {

        return this.patch.get(index);
    }

    public Double getFitness(int index) {

        return this.fitness.get(index);
    }

    public Integer getParent(Integer genotype) {

        int index = this.genotype.indexOf(genotype);
        return this.parent.get(index);
    }


    public List<Double> getCompleteBirths() {
        return this.birth;
    }
    public List<Double> getCompleteDeaths() {
        return this.death;
    }

    public List<Integer> getCompleteGenotypes() {
        return this.genotype;
    }



}
