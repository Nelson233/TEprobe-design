# -*- coding: utf-8 -*-
"""
Guideline for the installation is available  the website:https://docs.nupack.org/start/
Operation system：Mac/Linux operating systems or on the Linux subsystem of Windows 10/11

Users should update DNA_a1 and DNA_a2 to be the SNV-containing input RNA sequence and SNV-lacking input RNA sequence 
"""

from nupack import *


def CalcEnergy(DNA_a1, DNA_a2, CutIndex=0):
    
#Construction of TEprobes with toeholds
    Toehold_rec = 'GATAGGC'
    
    Toehold_block = 'CCCTATC'
    
    DNA_b_tmp = DNA_a1.replace('A','t').replace('T','a').replace('C', 'g').replace('G','c')
    
    DNA_b_tmp = DNA_b_tmp[::-1]
    
    DNA_b = Toehold_rec + DNA_b_tmp
    
    DNA_c = DNA_a1 + Toehold_block
    
    DNA_c = DNA_c[CutIndex:]

#Specify each strand involved in TEprobes and target sequences    
    a1 = Strand(DNA_a1, name='a1')
    
    a2 = Strand(DNA_a2, name='a2')
    
    b = Strand(DNA_b, name='b')
    
    c = Strand(DNA_c, name='c')
    
    d = Complex([a1, b], name='d')
    
    e = Complex([a2, b], name='e')
    
    f = Complex([b, c], name='f')

#Specify a physical model  
    model1 = Model(material='dna', ensemble='stacking', celsius=25, sodium=1.0, magnesium=0.0)
    
    set1 = ComplexSet(strands=[a1,b], complexes=SetSpec(max_size=2, include=[d]))
    
    set2 = ComplexSet(strands=[a2,b], complexes=SetSpec(max_size=2, include=[e]))
    
    set3 = ComplexSet(strands=[b,c], complexes=SetSpec(max_size=2, include=[f]))
     
#Calculate the free energy  
    complex_results1 = complex_analysis(complexes=set1, model=model1, compute=['mfe'])
    
    complex_results2 = complex_analysis(complexes=set2, model=model1, compute=['mfe'])
    
    complex_results3 = complex_analysis(complexes=set3, model=model1, compute=['mfe'])
    
    d_result = complex_results1[d]
    
    e_result = complex_results2[e]
    
    f_result = complex_results3[f]
    
    energy_d_result = d_result.mfe[0].energy

    energy_e_result = e_result.mfe[0].energy

    energy_f_result = f_result.mfe[0].energy

    delta_G_SNV_c = energy_d_result-(energy_f_result - 5.97)

    delta_G_SNV_l = energy_e_result-(energy_f_result - 5.97)

    return [delta_G_SNV_c, delta_G_SNV_l, DNA_b, DNA_c]


if __name__ == '__main__':
    
    DNA_a1 = 'TAAATACCATCCCCACGGCGATTTCGCAGTGTATG'

    DNA_a2 = 'TAAATACCATCCCCACGGCGATTCCGCAGTGTATG'


    print('DNA_a1: {}\n'.format(DNA_a1))
    
    print('DNA_a2: {}\n'.format(DNA_a2))
    
    
    CutIndex = 0
    
    [delta_G_SNV_c, delta_G_SNV_l, DNA_b, DNA_c] = CalcEnergy(DNA_a1, DNA_a2, CutIndex)
    
    
    print('\nCutIndex: {}\n'.format(CutIndex))
    print('DNA_c: {}\n'.format(DNA_c))

    print('delta_G_SNV_c: {}\n'.format(delta_G_SNV_c))
    print('delta_G_SNV_l: {}\n'.format(delta_G_SNV_l))
    
    if delta_G_SNV_c > 0:
        CutIndex = CutIndex + 1
        
        [delta_G_SNV_c, delta_G_SNV_l, DNA_b, DNA_c] = CalcEnergy(DNA_a1, DNA_a2, CutIndex)
        
        print('\n============ Iterations ============\n')
        
        print('\nCutIndex: {}\n'.format(CutIndex))
        print('DNA_c: {}\n'.format(DNA_c))
        print('delta_G_SNV_c: {}\n'.format(delta_G_SNV_c))
        print('delta_G_SNV_l: {}\n'.format(delta_G_SNV_l))
        
        while True:
            if delta_G_SNV_c > -2 and delta_G_SNV_c < 0 and delta_G_SNV_l > 0:
                break
            
            CutIndex = CutIndex + 1
            [delta_G_SNV_c, delta_G_SNV_l, DNA_b, DNA_c] = CalcEnergy(DNA_a1, DNA_a2, CutIndex)
            
            print('\nCutIndex: {}\n'.format(CutIndex))
            print('DNA_c: {}\n'.format(DNA_c))
            print('delta_G_SNV_c: {}\n'.format(delta_G_SNV_c))
            print('delta_G_SNV_l: {}\n'.format(delta_G_SNV_l))


        print('\n============ Output Results ============\n')

        for idx in range(10):
            print("NO.{} Block: {}".format(idx+1, DNA_c[idx:]))
        
        print('\nRec: {}\n'.format(DNA_b)) #Output TEprobe candidates

        
    else:
        print('ERROR: delta_G_SNV_c <= 0!')
