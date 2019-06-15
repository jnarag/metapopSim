import java.util.ArrayList;
import java.util.List;

/**
 * Created by jayna on 18/05/2018.
 */
public class infectionHistory {

    List<Integer> genotype = null;
    List<Double> birth = null;
    List<Integer> parent = null;
    List<List<Integer>> prevalence = null;
    List<Integer> patch = null;
    List<Double> fitness = null;
    List<Integer> parentOrigin = null;

    public infectionHistory () {

        genotype = new ArrayList<>();
        birth = new ArrayList<>();
        parent = new ArrayList<>();
        prevalence = new ArrayList<>();
        patch = new ArrayList<>();
        fitness = new ArrayList<>();
        parentOrigin = new ArrayList<>();

    }

    public void logGenotype(Integer genotype) {

        this.genotype.add(genotype);

    }

    public void logBirth(Double birth) {

        this.birth.add(birth);
    }

    public void logParent(Integer parent) {

        this.parent.add(parent);
    }

    public void logPatch(Integer patch) {

        this.patch.add(patch);
    }

    public void logParentOrigin(Integer patchOrBlood) {
        // 0 for patch; 1 for blood
        this.parentOrigin.add(patchOrBlood);
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

    public Integer getPatch(int index) {

        return this.patch.get(index);
    }

    public Double getFitness(int index) {

        return this.fitness.get(index);
    }



}
