
class Genome {

    constructor() {
      this.inputNeuronGenes = []; 
      this.outputNeuronGenes = [];
      this.neuronGenes = [];  // input, output, or hidden
      this.synapseGenes = [];
    //   this.numLayers = 2;  // input and output layer by default
    }
  
    // initializeGenes takes a genePool and adds all the synapse genes (edges) and associated neuron genes
    // (source and target nodes) from genePool to the current genome. It creates a genome with default topology.
    initializeGenes(genePool) {
        // this.numLayers = genePool.layers.length;
      for (let i = 0; i < genePool.synapseGenes.length; i++) {
        this.addSynapseGene(genePool.synapseGenes[i]);
        this.addNeuronGene(genePool.synapseGenes[i].from);
        this.addNeuronGene(genePool.synapseGenes[i].to);
      }
    }

    // crossover creates a new genome from the two received genomes.
    // If both parents have the same synapse gene, inherit it. If only of them has it, inherit it randomly.
    static crossover(genome1, genome2) {
        var newGenome = new Genome();
        var allSynapseGenes = genome1.combineSynapseGenesNoRepeat(genome2);
        // order edges from lower layers to higher
        allSynapseGenes.sort(function(a,b){return a.layer - b.layer}); // check if in correct order!!!!!*******
        
        //console.log(allEdges);

        for (let i = 0; i < allSynapseGenes.length; i++) {
            // If both genomes have synapse i, inherit it.
            if (genome1.containsSynapseGene(allSynapseGenes[i]) && genome2.containsSynapseGene(allSynapseGenes[i])) {
                newGenome.addNeuronGene(allSynapseGenes[i].from);
                newGenome.addNeuronGene(allSynapseGenes[i].to);
                newGenome.addSynapseGene(allSynapseGenes[i]);
            } else if (random(1) < 0.5) { // Inherit randomly otherwise (flip a coin).
                if (newGenome.containsNeuronGene(allSynapseGenes[i].from)) { // Check if origin neuron exists, so we can add a new outgoing edge from it.
                    newGenome.addNeuronGene(allSynapseGenes[i].to);
                    newGenome.addSynapseGene(allSynapseGenes[i]); 
                }
            }
        }

        return newGenome;
    }

    // crossoverBiased creates a new genome using common edges from parents.
    // It assumes genome1 as the fittest genome and gives it higher priority on inheritance.
    static crossoverBiased(genome1, genome2) {
        var newGenome = new Genome();
        var allSynapseGenes = genome1.combineSynapseGenesNoRepeat(genome2);
        // order edges from lower layers to higher
        allSynapseGenes.sort(function(a,b){return a.layer - b.layer}); // check if in correct order!!!!!*******

        // inherit common edges
        for (let i = 0; i < allSynapseGenes.length; i++) {
            // If both genomes have synapse i, inherit it.
            if (genome1.containsSynapseGene(allSynapseGenes[i]) && genome2.containsSynapseGene(allSynapseGenes[i])) {
                newGenome.addSynapseGene(allSynapseGenes[i]);
                newGenome.addNeuronGene(allSynapseGenes[i].from);
                newGenome.addNeuronGene(allSynapseGenes[i].to);
            } else if (genome1.containsSynapseGene(allSynapseGenes[i])) { // Prioritize synapse genes from genome1. Second pass to genes from parent1.
                if (random(1) < 0.7) { // threshold for passing this trait
                    if (newGenome.neuronExists(allSynapseGenes[i].from)) { // Check if origin neuron exists, so we can add a new outgoing edge from it.
                        newGenome.addNeuronGene(allSynapseGenes[i].to);
                        newGenome.addSynapseGene(allSynapseGenes[i]);
                    }
                }
            }
        }

        return newGenome;
    }

    // mutate
    mutate(newSynapseMutationRate, newNeuronMutationRate, genePool) {
        console.log("within mutation: synMutRate " + newSynapseMutationRate + ", neuMutRate " + newNeuronMutationRate);

        if (random(1) < newSynapseMutationRate) {
            // select randomly two nodes from my neuron genes
            let allNeurons = this.getAllNeuronGenes();
            // n1 and n2 neuron genes belong to my genome
            let n1 = random(allNeurons);
            let n2 = null;
            do {
                n2 = random(allNeurons);
            } while(n1.id == n2.id);

            //console.log(this.co);
            console.log("before mutation, to be joined (" + n1.id + "," + n2.id + ")");
            console.log(this);
            this.mutateAddSynapse(n1, n2, genePool);
            console.log("after mutation");
            console.log(this);
        }
        if (random(1) < newNeuronMutationRate) {
            let s1 = random(this.synapseGenes);

            console.log(this.containsSynapseGene(s1)); // should be true always
            console.log("before mutation, to be split id " + s1.id + " (" + s1.from.id + "," + s1.to.id + ")");
            console.log(this);
            this.mutateAddNeuron(s1, genePool);
            console.log("after mutation");
            console.log(this);
            console.log(this.containsSynapseGene(s1));
        }
    }

    // mutateAddSynapse looks for the (ng1, ng2) synapse in my genome and if it does not exist, it is created and
    // added to both genePool and my genome.
    // It receives NeuronGene objects ng1 and ng2 as parameters.
    mutateAddSynapse(ng1, ng2, genePool) {
        if (ng1.layer == ng2.layer) {
            return;
        } 
        if (ng1.layer > ng2.layer) {
            let temp = ng2;
            ng2 = ng1;
            ng1 = temp;
        } 
        // Check if synapse gene already exists in gene pool.
        let sg = genePool.getSynapseGene(ng1, ng2);
        // Create (ng1, ng2) synapse if it does not exist and add it to the gene pool.
        if (sg === null) { 
            let newSg = new SynapseGene(ng1, ng2);
            // New mutations are added to gene pool first
            genePool.addSynapseGene(newSg);
            this.addSynapseGene(newSg);
        } else {
            // Add synapse gene if it does not already exist in current genome.
            if (!this.containsSynapseGene(sg)) {
                this.addSynapseGene(sg);
            }
        }

        console.log("Successfully addded synapse mutation!");
    }

    // mutateAddNeuron receives a synapse gene sg and adds a new neuron in between sg.from and sg.to.
    // It adds the new edge to 
    mutateAddNeuron(sg, genePool) {
        var newLayerIndex = sg.from.layer + 1;
        // split edge
        var ng1 = sg.from;
        var ng2 = sg.to;

        // If edge crosses 1 layer, add new layer in between.
        if (sg.length === 1) {
            genePool.addLayer(newLayerIndex); // addLayer method updates the layer numbering for all nodes on layer > newLayerIndex
        }

        // Remove sg from my synapseGenes.
        this.removeSynapseGene(sg);
        // Add new node in between ng1 and ng2.
        // All neurons created at this point are of type hidden.
        var newNg = new NeuronGene(newLayerIndex, NODETYPE.HIDDEN);
        genePool.addNeuronGene(newNg);
        this.addNeuronGene(newNg);
        // create new edges
        var sg1 = new SynapseGene(ng1, newNg);
        var sg2 = new SynapseGene(newNg, ng2);
        genePool.addSynapseGene(sg1);
        genePool.addSynapseGene(sg2);
        this.addSynapseGene(sg1);
        this.addSynapseGene(sg2);

        console.log("Successfully added neuron mutation!!");
    }

    // It assumes each single neuron is only present in one of the subsets (does not check for repeats).
    getAllNeuronGenes() {
        var allNeuronGenes = [];
        for (let itm of this.inputNeuronGenes) {
            allNeuronGenes.push(itm);
        }
        for (let itm of this.outputNeuronGenes) {
            allNeuronGenes.push(itm);
        }
        for (let itm of this.neuronGenes) {
            allNeuronGenes.push(itm);
        }
        return allNeuronGenes;
    }
    
    addNeuronGene(ng) {
        if (!this.containsNeuronGene(ng)) { // avoids duplicates
            if (ng.type == NODETYPE.INPUT) {
                this.inputNeuronGenes.push(ng);
            } else if (ng.type == NODETYPE.OUTPUT) {
                this.outputNeuronGenes.push(ng);
            } else {
                this.neuronGenes.push(ng); // check only contains hidden type neurons!!!!!!!***********
            }
        }
    }
  
    addSynapseGene(sg) {
        if (!this.containsSynapseGene(sg)) { // avoids duplicates
            this.synapseGenes.push(sg);
        }
    }

    removeSynapseGene(sg) {
        for (let i = this.synapseGenes.length - 1; i >= 0; i--) {
            if (sg.id === this.synapseGenes[i].id) {
                this.synapseGenes.splice(i, 1);
            }
        }
    }

    containsSynapseGene(sg) {
        for (let i = 0; i < this.synapseGenes.length; i++) {
            if (sg.id === this.synapseGenes[i].id) {
                return true;
            }
        }
        return false;
    }

    containsNeuronGene(ng) {
        for (let i = 0; i < this.inputNeuronGenes.length; i++) {
            if (ng.id === this.inputNeuronGenes[i].id) {
                return true;
            }
        }
        for (let i = 0; i < this.outputNeuronGenes.length; i++) {
            if (ng.id === this.outputNeuronGenes[i].id) {
                return true;
            }
        }
        for (let i = 0; i < this.neuronGenes.length; i++) {
            if (ng.id === this.neuronGenes[i].id) {
                return true;
            }
        }
        return false;
    }

    combineSynapseGenesNoRepeat(genome) {
        var combined = [];
        combined = combined.concat(this.synapseGenes);

        for (let i = 0; i < genome.synapseGenes.length; i++) {
            let present = false;
            for (let j = 0; j < combined.length; j++) {
                if (genome.synapseGenes[i].id === combined[j].id) {
                    present = true;
                    break;
                }
            }
            if (!present) {
                combined.push(genome.synapseGenes[i]);
            }
        }
        return combined;
    }

    combineNeuronGenesNoRepeat(genome) {
        var combined = [];
        combined = combined.concat(this.getAllNeuronGenes());

        var gNeurons = genome.getAllNeuronGenes();
        for (let i = 0; i < gNeurons.length; i++) {
            if (!this.containsNeuronGene(gNeurons[i])) {
                combined.push(gNeurons[i]);
            }
        }
        return combined;
    }

    getSynapseGeneById(id) {
        for (let i = 0; i < this.synapseGenes.length; i++) {
          if (this.synapseGenes[i].id == id) {
            return s;
          }
        }
        return null;
    }

    getMaxSynapseLength() {
        var maxLength = 0;
        for (let s of this.synapseGenes) {
            if (s.length > maxLength) {
                maxLength = s.length;
            }
        }
        return maxLength;
    }
  
    copy() { // when is this used??????????
      var g = new Genome();
      for (let n of this.neuronGenes) {
        g.addNeuronGene(n.copy());
      }
      for (let c of this.synapseGenes) {
        g.addSynapseGene(c.copy());
      }
  
      return g;
    }
  
  }