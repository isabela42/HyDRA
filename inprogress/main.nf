// main.nf
workflow step01S1 {
    include { step01S1 } from './modules/step01S1.nf'

    step01S1()
}

workflow step03S1 {
    include { step03S1 } from './modules/step03S1.nf'

    step03S1()
}

