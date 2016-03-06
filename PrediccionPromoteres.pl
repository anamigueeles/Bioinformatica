#!/usr/bin/perl -w
# prog1.1
# Bruno Contreras-Moreira
#Citlali Gil Aguillon
#Jessica Danielly Medina Sánchez
#Anali Migueles Lozano
# Nearest Neighbor dG calculator

use strict;
sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }
# global variables
#variable donde guardamos primera mitad de la ventana para poder ver correccion simetrica
my $seqWin='';
#variable en la que se guarda segunda mitad de ventana traducida para ver correccion simetrica
my $seqWinR='';
#vector en el que se guardan las energias libres de todos los dinucleotidos 
my @dGn;
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my %NNparams    = (
        # SantaLucia J (1998) PNAS 95(4): 1460-1465.
        # [NaCl] 1M, 37C & pH=7
        # H(enthalpy): kcal/mol , S(entropy): cal/k.mol
        # stacking dinucleotides
        'AA/TT' , {'H',-7.9, 'S',-22.2},
        'AT/TA' , {'H',-7.2, 'S',-20.4},
        'TA/AT' , {'H',-7.2, 'S',-21.3},
        'CA/GT' , {'H',-8.5, 'S',-22.7},
        'GT/CA' , {'H',-8.4, 'S',-22.4},
        'CT/GA' , {'H',-7.8, 'S',-21.0},
        'GA/CT' , {'H',-8.2, 'S',-22.2},
        'CG/GC' , {'H',-10.6,'S',-27.2},
        'GC/CG' , {'H',-9.8, 'S',-24.4},
        'GG/CC' , {'H',-8.0, 'S',-19.9},
        # initiation costs
        'G'     , {'H', 0.1, 'S',-2.8 },
        'A'     , {'H', 2.3, 'S',4.1  },
        # symmetry correction
        'sym'   , {'H',   0, 'S',-1.4 } );
my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";
print "# parameters: Temperature=$T \c Window=$windowL\n\n";

open(SEQ, $infile) || die "# cannot open input $infile : $!\n ";
while(<SEQ>){
        if(/^(b\d{4}) \\ ([ATGC]+)/){
                my ($name,$seq) = ($1,$2);
                printf("sequence %s : ",$name);
						 #####BLOQUE PRIMERO#####

                ##MODIFICACION CODIGO: SUBRUTINA ACOPLADA CON CODIGO
                # calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
                # parameters: 1) DNA sequence string; 2) Celsius temperature
                # returns; 1) free energy scalar
                # uses global hash %NNparams

                #variables DNAstep: contiene dinucleotido,  nt: utilizado para correciones,  dG: se utiliza como intermediario en energia libre
                my ($DNAstep,$nt,$dG) = ('','',0);
                my @sequence = split(//,uc($seq));
                my $tK = 273.15 + $T;
                #vector en el que contendra el valor que corresponde a cada ventana
                my @dGn=();
                #creacion de VEctOR  dG, donde se guardaran para cada 15 sitios de la secuencia
                # add dG for overlapping dinculeotides

            	#loop que recorre secuencia para las posibles ventanas, por tanto lo hace 436 veces
                for (my $n=0; $n<(length($seq)-14); $n++){
            	           	#loop anidado que recorre las 15 posiciones siguiente formando la ventana
                            for(my $i=0; $i<$windowL-1; $i++){
                                $dG=0;		
                                $DNAstep='';
                                ####CAMBIO
                                    #ya que la subrutina afecta a la secuencia introducida, debemos usar variables temporales para no afectar la secuencia original
                                my $temporal=$sequence[$n+$i].$sequence[$n+$i+1];
                                $temporal= complement($temporal);
                                $DNAstep = $sequence[$n+$i].$sequence[$n+$i+1].'/'.$temporal;

                                if(!defined($NNparams{$DNAstep})){
                                    $DNAstep = reverse($DNAstep);
                                }

                                $dG = ((1000*$NNparams{$DNAstep}{'H'})-($tK*$NNparams{$DNAstep}{'S'}))/ 1000;
                                #rellenar vector
                                $dGn[$n] += $dG;

                                 #realizar bloque solo cuando se encuentre en la ultima vuelta
                                if($i==$windowL-2){

                                    # add correction for helix initiation
                                    $nt = $sequence[$n]; # first pair
                                    if(!defined($NNparams{$nt})){ 
                                        $nt = complement($nt) 
                                    }
                                    $dGn[$n] += ((1000*$NNparams{$nt}{'H'})-($tK*$NNparams{$nt}{'S'}))/ 1000;

                                    $nt = $sequence[$n+14]; # last pair
                                    if(!defined($NNparams{$nt})){ $
                                        nt = complement($nt) 
                                    }
                                    $dGn[$n] += ((1000*$NNparams{$nt}{'H'})-($tK*$NNparams{$nt}{'S'}))/ 1000;

                                    #please complete for symmetry correction [ AD ]
                                    #checar si la secuencia es palindromica  CAMBIOO
				                    for(my $a=0; $a<7; $a++){
                                        $seqWin .= $sequence[$n+$a];
                                        $seqWinR .= $sequence[$n+14-$a];
                                    }
                                    #sacar complementaria de la segunda parte de la ventana
                                    $seqWinR= complement($seqWinR);
                                    #si se cumple la igualdad debe sacarse correccion para ella
                                    if($seqWin eq $seqWinR){
                                        #se agrega a la ventana correspondiente
                                        $dGn[$n] += ((1000*$NNparams{'sym'}{'H'})-($tK*$NNparams{'sym'}{'S'}))/1000;
                                    }
                                }
                            }
            }
                                                        #######BLOQUE SEGUNDO######### 
                                                        # en donde se sacaran E1 E2 D#
#valores de cortes tomados del articulo de referencia
my $cutoff_1= 3.4;
my $cutoff_2= -15.99;
    ##CODIGO IMPLEMENTADO
    #calcular E1, E2, D
    #vector donde guardaremos posicion de ventana 
	my @dGSeq;
    #valores para E1
	my $E1=0; 
    #valores para E2
	my $E2=0;
    #valores para D
	my $D=0;
#recorrer todas las ventanas para calcular D
my $j=0;
#recorre todas las ventanas, con variable 'i'
for( my $i=0; $i<((scalar(@dGn))-199); $i++){
	#calcula E1
	for (my $n=0; $n<49; $n++){
		$E1+=$dGn[$n+$i];
	}#promedio 
    $E1/=50;
	#calcula E2
	for (my $n=0; $n<99; $n++){
		$E2+=$dGn[$n+99+$i];
	}#promedio
    $E2/=100;
	#diferencia de E1 y E2
	$D= $E1- $E2;
	#con el vector vamos a evaluar si respeta los puntos de corte
	if( $cutoff_2 > $E1 && $cutoff_1 > $D){
		#guarda posicion de valor relevante 
        ##de acuerdo con los criterios revisamos que se encuentre entre las posiciones -150 y 50 
		if(($i-350)>-150 && ($i-350)<50){
            #se guarda la posicion del posible promotor y se aumenta el contador 'j'
			$dGSeq[$j++]=$i-350;
            #se agregan 24 posiciones, porque los colindantes se toman como la misma señal (son 25, pero se aumenta otro en el 'for')
			$i+=24;		
		}
	}
    #reiniciar las variables para no tener problemas
	$E1=0;
	$E2=0;
	$D=0;
}
#imprimimos vector en el que se guardan los inicios de la secuencias predichas donde los valores son aceptables 
print join(", ", @dGSeq);
print "\n";
}
}close(SEQ);
print "\n";