const DISABLED = {
  newNeuron: 0,
  newSynapse: 0,
  randomWeight: 0,
  randomBias: 0,
  randomThreshold: 0
}

const SLOW = {
  newNeuron: 0.1,
  newSynapse: 0.1,
  randomWeight: 0.1,
  randomBias: 0.1,
  randomThreshold: 0.1
}

const MEDIUM = {
  newNeuron: 0.3,
  newSynapse: 0.3,
  randomWeight: 0.3,
  randomBias: 0.3,
  randomThreshold: 0.3
}

const FAST = {
  newNeuron: 0.5,
  newSynapse: 0.5,
  randomWeight: 0.5,
  randomBias: 0.5,
  randomThreshold: 0.5
}

class NeuronGene {

  constructor(layer, type) {
    this.id = 0;
    // this.threshold = t;
    // this.bias = b;
    this.layer = layer;
    this.type = type;
  }

  copy() { // When is this used????????
    return this;
  }
}

class SynapseGene {
  // constructor receives:
  // from: incoming NeuronGene object
  // to: outgoing NeuronGene object
  // ...
  constructor(from, to, enabled = true) {
    this.id = 0;
    this.from = from;
    this.to = to;
    this.enabled = enabled;
  }

  get length() {
    return this.to.layer - this.from.layer;
  }
 
  copy() { 
    return this;
  }

  disable() {
    this.enabled = false;
  }
}



