---
title: "Cluster Profiler Analysis for aTSC Mammary Glands"
author: "Dave Bridges"
date: "December 21, 2020"
output:
  html_document:
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
  pdf_document:
    highlight: tango
    keep_tex: yes
    number_sections: yes
    toc: yes
---



# Purpose

To use cluster profiler as part of DOSE to do gene set enrichments

# Raw Data

GSEA was run with folders put in this subfolder



# Gene Ontology


Table: Significant GO - BP Pathays

|           |ID         |Description                                                                                                                                      | pvalue| p.adjust|   NES|
|:----------|:----------|:------------------------------------------------------------------------------------------------------------------------------------------------|------:|--------:|-----:|
|GO:0042113 |GO:0042113 |B cell activation                                                                                                                                |  0.000|    0.000| -2.38|
|GO:0002250 |GO:0002250 |adaptive immune response                                                                                                                         |  0.000|    0.000| -2.26|
|GO:0030098 |GO:0030098 |lymphocyte differentiation                                                                                                                       |  0.000|    0.000| -2.10|
|GO:0050851 |GO:0050851 |antigen receptor-mediated signaling pathway                                                                                                      |  0.000|    0.000| -2.39|
|GO:0055002 |GO:0055002 |striated muscle cell development                                                                                                                 |  0.000|    0.000|  2.33|
|GO:0055001 |GO:0055001 |muscle cell development                                                                                                                          |  0.000|    0.000|  2.23|
|GO:1903131 |GO:1903131 |mononuclear cell differentiation                                                                                                                 |  0.000|    0.000| -1.96|
|GO:0051249 |GO:0051249 |regulation of lymphocyte activation                                                                                                              |  0.000|    0.000| -1.99|
|GO:0070661 |GO:0070661 |leukocyte proliferation                                                                                                                          |  0.000|    0.000| -2.10|
|GO:0002768 |GO:0002768 |immune response-regulating cell surface receptor signaling pathway                                                                               |  0.000|    0.000| -2.16|
|GO:0046651 |GO:0046651 |lymphocyte proliferation                                                                                                                         |  0.000|    0.000| -2.12|
|GO:0002443 |GO:0002443 |leukocyte mediated immunity                                                                                                                      |  0.000|    0.000| -2.05|
|GO:0032943 |GO:0032943 |mononuclear cell proliferation                                                                                                                   |  0.000|    0.000| -2.11|
|GO:0002366 |GO:0002366 |leukocyte activation involved in immune response                                                                                                 |  0.000|    0.000| -2.12|
|GO:0030239 |GO:0030239 |myofibril assembly                                                                                                                               |  0.000|    0.000|  2.39|
|GO:0002263 |GO:0002263 |cell activation involved in immune response                                                                                                      |  0.000|    0.000| -2.12|
|GO:0051146 |GO:0051146 |striated muscle cell differentiation                                                                                                             |  0.000|    0.000|  2.04|
|GO:0002429 |GO:0002429 |immune response-activating cell surface receptor signaling pathway                                                                               |  0.000|    0.000| -2.12|
|GO:0042692 |GO:0042692 |muscle cell differentiation                                                                                                                      |  0.000|    0.000|  1.93|
|GO:0002460 |GO:0002460 |adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains                        |  0.000|    0.000| -2.10|
|GO:0010927 |GO:0010927 |cellular component assembly involved in morphogenesis                                                                                            |  0.000|    0.000|  2.25|
|GO:0032989 |GO:0032989 |cellular anatomical entity morphogenesis                                                                                                         |  0.000|    0.000|  2.25|
|GO:0002377 |GO:0002377 |immunoglobulin production                                                                                                                        |  0.000|    0.000| -2.27|
|GO:0030183 |GO:0030183 |B cell differentiation                                                                                                                           |  0.000|    0.000| -2.27|
|GO:0002449 |GO:0002449 |lymphocyte mediated immunity                                                                                                                     |  0.000|    0.000| -2.09|
|GO:0002764 |GO:0002764 |immune response-regulating signaling pathway                                                                                                     |  0.000|    0.000| -1.85|
|GO:0050863 |GO:0050863 |regulation of T cell activation                                                                                                                  |  0.000|    0.000| -1.94|
|GO:0002285 |GO:0002285 |lymphocyte activation involved in immune response                                                                                                |  0.000|    0.000| -2.15|
|GO:0003012 |GO:0003012 |muscle system process                                                                                                                            |  0.000|    0.000|  1.87|
|GO:0046631 |GO:0046631 |alpha-beta T cell activation                                                                                                                     |  0.000|    0.000| -2.13|
|GO:0016064 |GO:0016064 |immunoglobulin mediated immune response                                                                                                          |  0.000|    0.000| -2.18|
|GO:0019724 |GO:0019724 |B cell mediated immunity                                                                                                                         |  0.000|    0.000| -2.17|
|GO:0007159 |GO:0007159 |leukocyte cell-cell adhesion                                                                                                                     |  0.000|    0.000| -1.88|
|GO:0002637 |GO:0002637 |regulation of immunoglobulin production                                                                                                          |  0.000|    0.000| -2.23|
|GO:0042100 |GO:0042100 |B cell proliferation                                                                                                                             |  0.000|    0.000| -2.19|
|GO:0051251 |GO:0051251 |positive regulation of lymphocyte activation                                                                                                     |  0.000|    0.000| -1.94|
|GO:0030217 |GO:0030217 |T cell differentiation                                                                                                                           |  0.000|    0.000| -1.91|
|GO:0031032 |GO:0031032 |actomyosin structure organization                                                                                                                |  0.000|    0.000|  1.98|
|GO:0070663 |GO:0070663 |regulation of leukocyte proliferation                                                                                                            |  0.000|    0.000| -1.97|
|GO:1903037 |GO:1903037 |regulation of leukocyte cell-cell adhesion                                                                                                       |  0.000|    0.000| -1.88|
|GO:0055013 |GO:0055013 |cardiac muscle cell development                                                                                                                  |  0.000|    0.000|  2.15|
|GO:0050853 |GO:0050853 |B cell receptor signaling pathway                                                                                                                |  0.000|    0.000| -2.25|
|GO:0050864 |GO:0050864 |regulation of B cell activation                                                                                                                  |  0.000|    0.000| -2.11|
|GO:0045214 |GO:0045214 |sarcomere organization                                                                                                                           |  0.000|    0.000|  2.22|
|GO:0055006 |GO:0055006 |cardiac cell development                                                                                                                         |  0.000|    0.000|  2.15|
|GO:0002757 |GO:0002757 |immune response-activating signaling pathway                                                                                                     |  0.000|    0.000| -1.79|
|GO:0014706 |GO:0014706 |striated muscle tissue development                                                                                                               |  0.000|    0.000|  1.76|
|GO:0050852 |GO:0050852 |T cell receptor signaling pathway                                                                                                                |  0.000|    0.000| -2.10|
|GO:0050670 |GO:0050670 |regulation of lymphocyte proliferation                                                                                                           |  0.000|    0.000| -1.95|
|GO:1903039 |GO:1903039 |positive regulation of leukocyte cell-cell adhesion                                                                                              |  0.000|    0.000| -1.93|
|GO:0002697 |GO:0002697 |regulation of immune effector process                                                                                                            |  0.000|    0.000| -1.82|
|GO:0002312 |GO:0002312 |B cell activation involved in immune response                                                                                                    |  0.000|    0.000| -2.15|
|GO:0032944 |GO:0032944 |regulation of mononuclear cell proliferation                                                                                                     |  0.000|    0.000| -1.95|
|GO:0035710 |GO:0035710 |CD4-positive, alpha-beta T cell activation                                                                                                       |  0.000|    0.000| -2.07|
|GO:0002562 |GO:0002562 |somatic diversification of immune receptors via germline recombination within a single locus                                                     |  0.000|    0.000| -2.18|
|GO:0016444 |GO:0016444 |somatic cell DNA recombination                                                                                                                   |  0.000|    0.000| -2.18|
|GO:0055007 |GO:0055007 |cardiac muscle cell differentiation                                                                                                              |  0.000|    0.000|  2.04|
|GO:0050870 |GO:0050870 |positive regulation of T cell activation                                                                                                         |  0.000|    0.000| -1.92|
|GO:0002440 |GO:0002440 |production of molecular mediator of immune response                                                                                              |  0.000|    0.000| -1.88|
|GO:0002253 |GO:0002253 |activation of immune response                                                                                                                    |  0.000|    0.000| -1.69|
|GO:0060537 |GO:0060537 |muscle tissue development                                                                                                                        |  0.000|    0.000|  1.70|
|GO:0002381 |GO:0002381 |immunoglobulin production involved in immunoglobulin-mediated immune response                                                                    |  0.000|    0.000| -2.17|
|GO:0002703 |GO:0002703 |regulation of leukocyte mediated immunity                                                                                                        |  0.000|    0.000| -1.88|
|GO:0070665 |GO:0070665 |positive regulation of leukocyte proliferation                                                                                                   |  0.000|    0.000| -1.98|
|GO:0050867 |GO:0050867 |positive regulation of cell activation                                                                                                           |  0.000|    0.000| -1.76|
|GO:0043500 |GO:0043500 |muscle adaptation                                                                                                                                |  0.000|    0.000|  2.04|
|GO:0002696 |GO:0002696 |positive regulation of leukocyte activation                                                                                                      |  0.000|    0.000| -1.79|
|GO:0002200 |GO:0002200 |somatic diversification of immune receptors                                                                                                      |  0.000|    0.000| -2.12|
|GO:0050854 |GO:0050854 |regulation of antigen receptor-mediated signaling pathway                                                                                        |  0.000|    0.000| -2.08|
|GO:0022409 |GO:0022409 |positive regulation of cell-cell adhesion                                                                                                        |  0.000|    0.000| -1.79|
|GO:0002520 |GO:0002520 |immune system development                                                                                                                        |  0.000|    0.000| -1.87|
|GO:0042098 |GO:0042098 |T cell proliferation                                                                                                                             |  0.000|    0.000| -1.86|
|GO:0032946 |GO:0032946 |positive regulation of mononuclear cell proliferation                                                                                            |  0.000|    0.000| -1.97|
|GO:0035051 |GO:0035051 |cardiocyte differentiation                                                                                                                       |  0.000|    0.000|  1.92|
|GO:0022407 |GO:0022407 |regulation of cell-cell adhesion                                                                                                                 |  0.000|    0.000| -1.66|
|GO:0002204 |GO:0002204 |somatic recombination of immunoglobulin genes involved in immune response                                                                        |  0.000|    0.000| -2.08|
|GO:0002208 |GO:0002208 |somatic diversification of immunoglobulins involved in immune response                                                                           |  0.000|    0.000| -2.08|
|GO:0045190 |GO:0045190 |isotype switching                                                                                                                                |  0.000|    0.000| -2.08|
|GO:0046632 |GO:0046632 |alpha-beta T cell differentiation                                                                                                                |  0.000|    0.000| -1.96|
|GO:0050671 |GO:0050671 |positive regulation of lymphocyte proliferation                                                                                                  |  0.000|    0.000| -1.98|
|GO:0045058 |GO:0045058 |T cell selection                                                                                                                                 |  0.000|    0.000| -2.06|
|GO:0032673 |GO:0032673 |regulation of interleukin-4 production                                                                                                           |  0.000|    0.000| -2.09|
|GO:0043501 |GO:0043501 |skeletal muscle adaptation                                                                                                                       |  0.000|    0.000|  2.06|
|GO:0016447 |GO:0016447 |somatic recombination of immunoglobulin gene segments                                                                                            |  0.000|    0.000| -2.06|
|GO:0002819 |GO:0002819 |regulation of adaptive immune response                                                                                                           |  0.000|    0.000| -1.85|
|GO:0006936 |GO:0006936 |muscle contraction                                                                                                                               |  0.000|    0.000|  1.73|
|GO:0043367 |GO:0043367 |CD4-positive, alpha-beta T cell differentiation                                                                                                  |  0.000|    0.000| -1.98|
|GO:0042129 |GO:0042129 |regulation of T cell proliferation                                                                                                               |  0.000|    0.000| -1.85|
|GO:0006941 |GO:0006941 |striated muscle contraction                                                                                                                      |  0.000|    0.000|  1.89|
|GO:0007517 |GO:0007517 |muscle organ development                                                                                                                         |  0.000|    0.000|  1.66|
|GO:0050871 |GO:0050871 |positive regulation of B cell activation                                                                                                         |  0.000|    0.000| -1.97|
|GO:0002821 |GO:0002821 |positive regulation of adaptive immune response                                                                                                  |  0.000|    0.000| -1.90|
|GO:0002683 |GO:0002683 |negative regulation of immune system process                                                                                                     |  0.000|    0.000| -1.58|
|GO:0002712 |GO:0002712 |regulation of B cell mediated immunity                                                                                                           |  0.000|    0.000| -2.01|
|GO:0002889 |GO:0002889 |regulation of immunoglobulin mediated immune response                                                                                            |  0.000|    0.000| -2.01|
|GO:0016445 |GO:0016445 |somatic diversification of immunoglobulins                                                                                                       |  0.000|    0.000| -2.01|
|GO:0042102 |GO:0042102 |positive regulation of T cell proliferation                                                                                                      |  0.000|    0.000| -1.92|
|GO:0002706 |GO:0002706 |regulation of lymphocyte mediated immunity                                                                                                       |  0.000|    0.000| -1.84|
|GO:0032633 |GO:0032633 |interleukin-4 production                                                                                                                         |  0.000|    0.000| -2.09|
|GO:0048738 |GO:0048738 |cardiac muscle tissue development                                                                                                                |  0.000|    0.000|  1.74|
|GO:0002639 |GO:0002639 |positive regulation of immunoglobulin production                                                                                                 |  0.000|    0.000| -2.03|
|GO:0002699 |GO:0002699 |positive regulation of immune effector process                                                                                                   |  0.000|    0.000| -1.75|
|GO:0002824 |GO:0002824 |positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains |  0.000|    0.000| -1.89|
|GO:0043368 |GO:0043368 |positive T cell selection                                                                                                                        |  0.000|    0.000| -2.03|
|GO:0002286 |GO:0002286 |T cell activation involved in immune response                                                                                                    |  0.000|    0.000| -1.87|
|GO:0002822 |GO:0002822 |regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains          |  0.000|    0.000| -1.79|
|GO:0140694 |GO:0140694 |non-membrane-bounded organelle assembly                                                                                                          |  0.000|    0.000|  1.59|
|GO:0002335 |GO:0002335 |mature B cell differentiation                                                                                                                    |  0.000|    0.000| -2.02|
|GO:0006310 |GO:0006310 |DNA recombination                                                                                                                                |  0.000|    0.000| -1.65|
|GO:0055003 |GO:0055003 |cardiac myofibril assembly                                                                                                                       |  0.000|    0.000|  1.98|
|GO:0001819 |GO:0001819 |positive regulation of cytokine production                                                                                                       |  0.000|    0.000| -1.55|
|GO:0002287 |GO:0002287 |alpha-beta T cell activation involved in immune response                                                                                         |  0.000|    0.001| -1.90|
|GO:0032753 |GO:0032753 |positive regulation of interleukin-4 production                                                                                                  |  0.000|    0.001| -2.01|
|GO:0002700 |GO:0002700 |regulation of production of molecular mediator of immune response                                                                                |  0.000|    0.001| -1.76|
|GO:0030888 |GO:0030888 |regulation of B cell proliferation                                                                                                               |  0.000|    0.001| -1.96|
|GO:0002456 |GO:0002456 |T cell mediated immunity                                                                                                                         |  0.000|    0.001| -1.84|
|GO:0045066 |GO:0045066 |regulatory T cell differentiation                                                                                                                |  0.000|    0.001| -2.02|
|GO:0002274 |GO:0002274 |myeloid leukocyte activation                                                                                                                     |  0.000|    0.001| -1.68|
|GO:0002638 |GO:0002638 |negative regulation of immunoglobulin production                                                                                                 |  0.000|    0.001| -1.92|
|GO:0002293 |GO:0002293 |alpha-beta T cell differentiation involved in immune response                                                                                    |  0.000|    0.001| -1.89|
|GO:0045619 |GO:0045619 |regulation of lymphocyte differentiation                                                                                                         |  0.000|    0.001| -1.70|
|GO:0043502 |GO:0043502 |regulation of muscle adaptation                                                                                                                  |  0.000|    0.001|  1.81|
|GO:0003009 |GO:0003009 |skeletal muscle contraction                                                                                                                      |  0.000|    0.001|  1.94|
|GO:0002292 |GO:0002292 |T cell differentiation involved in immune response                                                                                               |  0.000|    0.001| -1.86|
|GO:0002294 |GO:0002294 |CD4-positive, alpha-beta T cell differentiation involved in immune response                                                                      |  0.000|    0.001| -1.87|
|GO:0045785 |GO:0045785 |positive regulation of cell adhesion                                                                                                             |  0.000|    0.001| -1.50|
|GO:0050881 |GO:0050881 |musculoskeletal movement                                                                                                                         |  0.000|    0.002|  1.89|
|GO:1903307 |GO:1903307 |positive regulation of regulated secretory pathway                                                                                               |  0.000|    0.002| -1.88|
|GO:1902107 |GO:1902107 |positive regulation of leukocyte differentiation                                                                                                 |  0.000|    0.002| -1.67|
|GO:1903708 |GO:1903708 |positive regulation of hemopoiesis                                                                                                               |  0.000|    0.002| -1.67|
|GO:0045621 |GO:0045621 |positive regulation of lymphocyte differentiation                                                                                                |  0.000|    0.002| -1.78|
|GO:0045055 |GO:0045055 |regulated exocytosis                                                                                                                             |  0.000|    0.002| -1.63|
|GO:0042093 |GO:0042093 |T-helper cell differentiation                                                                                                                    |  0.000|    0.002| -1.85|
|GO:0046634 |GO:0046634 |regulation of alpha-beta T cell activation                                                                                                       |  0.000|    0.003| -1.77|
|GO:0048535 |GO:0048535 |lymph node development                                                                                                                           |  0.000|    0.003| -1.97|
|GO:0002705 |GO:0002705 |positive regulation of leukocyte mediated immunity                                                                                               |  0.000|    0.003| -1.73|
|GO:0010720 |GO:0010720 |positive regulation of cell development                                                                                                          |  0.000|    0.003| -1.45|
|GO:0002260 |GO:0002260 |lymphocyte homeostasis                                                                                                                           |  0.000|    0.003| -1.79|
|GO:0002708 |GO:0002708 |positive regulation of lymphocyte mediated immunity                                                                                              |  0.000|    0.003| -1.78|
|GO:2000514 |GO:2000514 |regulation of CD4-positive, alpha-beta T cell activation                                                                                         |  0.000|    0.003| -1.80|
|GO:0032609 |GO:0032609 |type II interferon production                                                                                                                    |  0.000|    0.004| -1.71|
|GO:0002714 |GO:0002714 |positive regulation of B cell mediated immunity                                                                                                  |  0.000|    0.004| -1.93|
|GO:0002891 |GO:0002891 |positive regulation of immunoglobulin mediated immune response                                                                                   |  0.000|    0.004| -1.93|
|GO:0050879 |GO:0050879 |multicellular organismal movement                                                                                                                |  0.000|    0.004|  1.87|
|GO:0043299 |GO:0043299 |leukocyte degranulation                                                                                                                          |  0.000|    0.004| -1.79|
|GO:0002275 |GO:0002275 |myeloid cell activation involved in immune response                                                                                              |  0.000|    0.004| -1.78|
|GO:0045580 |GO:0045580 |regulation of T cell differentiation                                                                                                             |  0.000|    0.004| -1.67|
|GO:0048291 |GO:0048291 |isotype switching to IgG isotypes                                                                                                                |  0.000|    0.004| -1.90|
|GO:0032660 |GO:0032660 |regulation of interleukin-17 production                                                                                                          |  0.000|    0.005| -1.86|
|GO:2000516 |GO:2000516 |positive regulation of CD4-positive, alpha-beta T cell activation                                                                                |  0.000|    0.005| -1.94|
|GO:0045191 |GO:0045191 |regulation of isotype switching                                                                                                                  |  0.000|    0.005| -1.92|
|GO:0071467 |GO:0071467 |cellular response to pH                                                                                                                          |  0.000|    0.005| -1.91|
|GO:0072678 |GO:0072678 |T cell migration                                                                                                                                 |  0.000|    0.005| -1.80|
|GO:0002361 |GO:0002361 |CD4-positive, CD25-positive, alpha-beta regulatory T cell differentiation                                                                        |  0.000|    0.005| -1.83|
|GO:1903706 |GO:1903706 |regulation of hemopoiesis                                                                                                                        |  0.000|    0.005| -1.49|
|GO:0045582 |GO:0045582 |positive regulation of T cell differentiation                                                                                                    |  0.000|    0.006| -1.72|
|GO:0002702 |GO:0002702 |positive regulation of production of molecular mediator of immune response                                                                       |  0.000|    0.006| -1.69|
|GO:0048302 |GO:0048302 |regulation of isotype switching to IgG isotypes                                                                                                  |  0.000|    0.006| -1.89|
|GO:0060538 |GO:0060538 |skeletal muscle organ development                                                                                                                |  0.000|    0.006|  1.60|
|GO:0090257 |GO:0090257 |regulation of muscle system process                                                                                                              |  0.000|    0.006|  1.57|
|GO:0002438 |GO:0002438 |acute inflammatory response to antigenic stimulus                                                                                                |  0.000|    0.008| -1.89|
|GO:0007259 |GO:0007259 |cell surface receptor signaling pathway via JAK-STAT                                                                                             |  0.000|    0.008| -1.68|
|GO:0014888 |GO:0014888 |striated muscle adaptation                                                                                                                       |  0.000|    0.008|  1.82|
|GO:0043320 |GO:0043320 |natural killer cell degranulation                                                                                                                |  0.000|    0.008| -1.83|
|GO:0046635 |GO:0046635 |positive regulation of alpha-beta T cell activation                                                                                              |  0.000|    0.008| -1.81|
|GO:0003007 |GO:0003007 |heart morphogenesis                                                                                                                              |  0.000|    0.008|  1.54|
|GO:0031294 |GO:0031294 |lymphocyte costimulation                                                                                                                         |  0.000|    0.008| -1.89|
|GO:0032729 |GO:0032729 |positive regulation of type II interferon production                                                                                             |  0.000|    0.008| -1.74|
|GO:0001776 |GO:0001776 |leukocyte homeostasis                                                                                                                            |  0.000|    0.008| -1.68|
|GO:0048266 |GO:0048266 |behavioral response to pain                                                                                                                      |  0.000|    0.008| -1.87|
|GO:0001906 |GO:0001906 |cell killing                                                                                                                                     |  0.000|    0.009| -1.63|
|GO:0014874 |GO:0014874 |response to stimulus involved in regulation of muscle adaptation                                                                                 |  0.000|    0.009|  1.87|
|GO:0030890 |GO:0030890 |positive regulation of B cell proliferation                                                                                                      |  0.000|    0.009| -1.85|
|GO:0032655 |GO:0032655 |regulation of interleukin-12 production                                                                                                          |  0.000|    0.009| -1.77|
|GO:0045348 |GO:0045348 |positive regulation of MHC class II biosynthetic process                                                                                         |  0.000|    0.010| -1.84|
|GO:0006909 |GO:0006909 |phagocytosis                                                                                                                                     |  0.000|    0.010| -1.57|
|GO:0007389 |GO:0007389 |pattern specification process                                                                                                                    |  0.000|    0.010|  1.44|
|GO:0032615 |GO:0032615 |interleukin-12 production                                                                                                                        |  0.000|    0.011| -1.77|
|GO:1903305 |GO:1903305 |regulation of regulated secretory pathway                                                                                                        |  0.000|    0.011| -1.63|
|GO:0097696 |GO:0097696 |cell surface receptor signaling pathway via STAT                                                                                                 |  0.000|    0.012| -1.62|
|GO:0060047 |GO:0060047 |heart contraction                                                                                                                                |  0.000|    0.012|  1.54|
|GO:0045589 |GO:0045589 |regulation of regulatory T cell differentiation                                                                                                  |  0.000|    0.012| -1.85|
|GO:0050855 |GO:0050855 |regulation of B cell receptor signaling pathway                                                                                                  |  0.000|    0.012| -1.84|
|GO:1902105 |GO:1902105 |regulation of leukocyte differentiation                                                                                                          |  0.000|    0.012| -1.49|
|GO:0007519 |GO:0007519 |skeletal muscle tissue development                                                                                                               |  0.000|    0.012|  1.60|
|GO:0043300 |GO:0043300 |regulation of leukocyte degranulation                                                                                                            |  0.000|    0.012| -1.80|
|GO:0043374 |GO:0043374 |CD8-positive, alpha-beta T cell differentiation                                                                                                  |  0.000|    0.013| -1.85|
|GO:0002861 |GO:0002861 |regulation of inflammatory response to antigenic stimulus                                                                                        |  0.000|    0.014| -1.86|
|GO:0031295 |GO:0031295 |T cell costimulation                                                                                                                             |  0.000|    0.014| -1.84|
|GO:0045346 |GO:0045346 |regulation of MHC class II biosynthetic process                                                                                                  |  0.000|    0.014| -1.84|
|GO:0009988 |GO:0009988 |cell-cell recognition                                                                                                                            |  0.000|    0.014| -1.79|
|GO:0050848 |GO:0050848 |regulation of calcium-mediated signaling                                                                                                         |  0.000|    0.014| -1.70|
|GO:0050856 |GO:0050856 |regulation of T cell receptor signaling pathway                                                                                                  |  0.000|    0.016| -1.77|
|GO:0140448 |GO:0140448 |signaling receptor ligand precursor processing                                                                                                   |  0.001|    0.016|  1.85|
|GO:0060443 |GO:0060443 |mammary gland morphogenesis                                                                                                                      |  0.001|    0.016|  1.79|
|GO:0019221 |GO:0019221 |cytokine-mediated signaling pathway                                                                                                              |  0.001|    0.016| -1.40|
|GO:0046633 |GO:0046633 |alpha-beta T cell proliferation                                                                                                                  |  0.001|    0.017| -1.77|
|GO:0001782 |GO:0001782 |B cell homeostasis                                                                                                                               |  0.001|    0.018| -1.77|
|GO:0060972 |GO:0060972 |left/right pattern formation                                                                                                                     |  0.001|    0.019|  1.67|
|GO:0072676 |GO:0072676 |lymphocyte migration                                                                                                                             |  0.001|    0.019| -1.68|
|GO:0032649 |GO:0032649 |regulation of type II interferon production                                                                                                      |  0.001|    0.019| -1.62|
|GO:0007600 |GO:0007600 |sensory perception                                                                                                                               |  0.001|    0.019| -1.40|
|GO:0048534 |GO:0048534 |hematopoietic or lymphoid organ development                                                                                                      |  0.001|    0.020| -1.64|
|GO:0002701 |GO:0002701 |negative regulation of production of molecular mediator of immune response                                                                       |  0.001|    0.020| -1.75|
|GO:0045059 |GO:0045059 |positive thymic T cell selection                                                                                                                 |  0.001|    0.022| -1.82|
|GO:0043372 |GO:0043372 |positive regulation of CD4-positive, alpha-beta T cell differentiation                                                                           |  0.001|    0.022| -1.80|
|GO:0007368 |GO:0007368 |determination of left/right symmetry                                                                                                             |  0.001|    0.022|  1.65|
|GO:0051153 |GO:0051153 |regulation of striated muscle cell differentiation                                                                                               |  0.001|    0.022|  1.65|
|GO:0009615 |GO:0009615 |response to virus                                                                                                                                |  0.001|    0.023| -1.45|
|GO:0045342 |GO:0045342 |MHC class II biosynthetic process                                                                                                                |  0.001|    0.024| -1.79|
|GO:0046637 |GO:0046637 |regulation of alpha-beta T cell differentiation                                                                                                  |  0.001|    0.024| -1.70|
|GO:0003015 |GO:0003015 |heart process                                                                                                                                    |  0.001|    0.025|  1.49|
|GO:0010324 |GO:0010324 |membrane invagination                                                                                                                            |  0.001|    0.025| -1.68|
|GO:0022612 |GO:0022612 |gland morphogenesis                                                                                                                              |  0.001|    0.025|  1.59|
|GO:0042509 |GO:0042509 |regulation of tyrosine phosphorylation of STAT protein                                                                                           |  0.001|    0.025| -1.69|
|GO:0034244 |GO:0034244 |negative regulation of transcription elongation by RNA polymerase II                                                                             |  0.001|    0.027|  1.82|
|GO:0002369 |GO:0002369 |T cell cytokine production                                                                                                                       |  0.001|    0.027| -1.78|
|GO:0033275 |GO:0033275 |actin-myosin filament sliding                                                                                                                    |  0.001|    0.027|  1.75|
|GO:0042330 |GO:0042330 |taxis                                                                                                                                            |  0.001|    0.027| -1.38|
|GO:0002313 |GO:0002313 |mature B cell differentiation involved in immune response                                                                                        |  0.001|    0.027| -1.79|
|GO:0000018 |GO:0000018 |regulation of DNA recombination                                                                                                                  |  0.001|    0.027| -1.59|
|GO:0008037 |GO:0008037 |cell recognition                                                                                                                                 |  0.001|    0.028| -1.58|
|GO:0046638 |GO:0046638 |positive regulation of alpha-beta T cell differentiation                                                                                         |  0.001|    0.028| -1.73|
|GO:0001708 |GO:0001708 |cell fate specification                                                                                                                          |  0.001|    0.028|  1.64|
|GO:0002920 |GO:0002920 |regulation of humoral immune response                                                                                                            |  0.001|    0.029| -1.78|
|GO:0050900 |GO:0050900 |leukocyte migration                                                                                                                              |  0.001|    0.029| -1.42|
|GO:0000272 |GO:0000272 |polysaccharide catabolic process                                                                                                                 |  0.001|    0.030|  1.80|
|GO:0051147 |GO:0051147 |regulation of muscle cell differentiation                                                                                                        |  0.001|    0.030|  1.57|
|GO:0071216 |GO:0071216 |cellular response to biotic stimulus                                                                                                             |  0.001|    0.031| -1.50|
|GO:0036230 |GO:0036230 |granulocyte activation                                                                                                                           |  0.001|    0.031| -1.76|
|GO:0050862 |GO:0050862 |positive regulation of T cell receptor signaling pathway                                                                                         |  0.001|    0.031| -1.82|
|GO:0019722 |GO:0019722 |calcium-mediated signaling                                                                                                                       |  0.001|    0.032| -1.50|
|GO:0016052 |GO:0016052 |carbohydrate catabolic process                                                                                                                   |  0.001|    0.032|  1.56|
|GO:0060742 |GO:0060742 |epithelial cell differentiation involved in prostate gland development                                                                           |  0.001|    0.033|  1.81|
|GO:0002444 |GO:0002444 |myeloid leukocyte mediated immunity                                                                                                              |  0.001|    0.033| -1.58|
|GO:0060603 |GO:0060603 |mammary gland duct morphogenesis                                                                                                                 |  0.001|    0.034|  1.73|
|GO:0061180 |GO:0061180 |mammary gland epithelium development                                                                                                             |  0.001|    0.034|  1.63|
|GO:0032620 |GO:0032620 |interleukin-17 production                                                                                                                        |  0.001|    0.034| -1.76|
|GO:0032640 |GO:0032640 |tumor necrosis factor production                                                                                                                 |  0.001|    0.036| -1.53|
|GO:0071706 |GO:0071706 |tumor necrosis factor superfamily cytokine production                                                                                            |  0.001|    0.037| -1.55|
|GO:0050906 |GO:0050906 |detection of stimulus involved in sensory perception                                                                                             |  0.001|    0.037| -1.59|
|GO:0009251 |GO:0009251 |glucan catabolic process                                                                                                                         |  0.001|    0.037|  1.75|
|GO:1903555 |GO:1903555 |regulation of tumor necrosis factor superfamily cytokine production                                                                              |  0.001|    0.038| -1.53|
|GO:0038093 |GO:0038093 |Fc receptor signaling pathway                                                                                                                    |  0.002|    0.039| -1.78|
|GO:0014733 |GO:0014733 |regulation of skeletal muscle adaptation                                                                                                         |  0.002|    0.039|  1.78|
|GO:0043370 |GO:0043370 |regulation of CD4-positive, alpha-beta T cell differentiation                                                                                    |  0.002|    0.039| -1.72|
|GO:0043302 |GO:0043302 |positive regulation of leukocyte degranulation                                                                                                   |  0.002|    0.040| -1.79|
|GO:0048665 |GO:0048665 |neuron fate specification                                                                                                                        |  0.002|    0.040|  1.75|
|GO:0002709 |GO:0002709 |regulation of T cell mediated immunity                                                                                                           |  0.002|    0.040| -1.61|
|GO:0032785 |GO:0032785 |negative regulation of DNA-templated transcription, elongation                                                                                   |  0.002|    0.041|  1.78|
|GO:0036336 |GO:0036336 |dendritic cell migration                                                                                                                         |  0.002|    0.041| -1.76|
|GO:0050857 |GO:0050857 |positive regulation of antigen receptor-mediated signaling pathway                                                                               |  0.002|    0.041| -1.76|
|GO:0050850 |GO:0050850 |positive regulation of calcium-mediated signaling                                                                                                |  0.002|    0.041| -1.65|
|GO:0043588 |GO:0043588 |skin development                                                                                                                                 |  0.002|    0.042|  1.43|
|GO:0002685 |GO:0002685 |regulation of leukocyte migration                                                                                                                |  0.002|    0.042| -1.49|
|GO:0048268 |GO:0048268 |clathrin coat assembly                                                                                                                           |  0.002|    0.042| -1.78|
|GO:0021513 |GO:0021513 |spinal cord dorsal/ventral patterning                                                                                                            |  0.002|    0.043|  1.75|
|GO:0048304 |GO:0048304 |positive regulation of isotype switching to IgG isotypes                                                                                         |  0.002|    0.043| -1.75|
|GO:0045061 |GO:0045061 |thymic T cell selection                                                                                                                          |  0.002|    0.044| -1.76|
|GO:0033555 |GO:0033555 |multicellular organismal response to stress                                                                                                      |  0.002|    0.045| -1.55|
|GO:1990868 |GO:1990868 |response to chemokine                                                                                                                            |  0.002|    0.045| -1.63|
|GO:1990869 |GO:1990869 |cellular response to chemokine                                                                                                                   |  0.002|    0.045| -1.63|
|GO:0002698 |GO:0002698 |negative regulation of immune effector process                                                                                                   |  0.002|    0.045| -1.60|
|GO:0006887 |GO:0006887 |exocytosis                                                                                                                                       |  0.002|    0.045| -1.39|
|GO:0007260 |GO:0007260 |tyrosine phosphorylation of STAT protein                                                                                                         |  0.002|    0.046| -1.64|
|GO:0071219 |GO:0071219 |cellular response to molecule of bacterial origin                                                                                                |  0.002|    0.046| -1.48|
|GO:0051321 |GO:0051321 |meiotic cell cycle                                                                                                                               |  0.002|    0.047| -1.44|
|GO:0070098 |GO:0070098 |chemokine-mediated signaling pathway                                                                                                             |  0.002|    0.047| -1.66|
|GO:0018212 |GO:0018212 |peptidyl-tyrosine modification                                                                                                                   |  0.002|    0.047| -1.44|
|GO:0045830 |GO:0045830 |positive regulation of isotype switching                                                                                                         |  0.002|    0.048| -1.75|
|GO:0033151 |GO:0033151 |V(D)J recombination                                                                                                                              |  0.002|    0.048| -1.75|
|GO:0002886 |GO:0002886 |regulation of myeloid leukocyte mediated immunity                                                                                                |  0.002|    0.048| -1.67|
|GO:0039692 |GO:0039692 |single stranded viral RNA replication via double stranded DNA intermediate                                                                       |  0.002|    0.049| -1.74|
|GO:0018108 |GO:0018108 |peptidyl-tyrosine phosphorylation                                                                                                                |  0.002|    0.049| -1.44|

![GO-BP Enrichment Plot](figures-noura/go-bp-1.png)![GO-BP Enrichment Plot](figures-noura/go-bp-2.png)

Table: Summary pf GO-BP pathways

|Significant |Direction     |   n|
|:-----------|:-------------|---:|
|Significant |Downregulated | 216|
|Significant |Upregulated   |  58|

# Wikipathways

From DisGeNET




# Cell Markers

```{ cell-markers, fig.cap="Cell markers enrichment plot"}
cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt')

cells <- cell_marker_data %>%
    dplyr::select(cellName, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
    tidyr::unnest()
cm <- GSEA(genes, TERM2GENE = cells) 

cm %>% select(ID, Description, pvalue,p.adjust, NES ) %>% kable(caption="Significant cell types")

cm.file <- 'GSEA Results - Cell Marker.csv'
cm %>%
  as.data.frame %>%
  write_csv(file=cm.file)

cm <- pairwise_termsim(cm, method="JC")
emapplot(cm, color='NES')
upsetplot(cm)
ridgeplot(cm, fill="NES")
```

# Transcription Factors



# Session Information


``` r
sessionInfo()
```

```
## R version 4.4.2 (2024-10-31)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.1.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Detroit
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] msigdbr_7.5.1          fgsea_1.32.0           ggnewscale_0.5.0      
##  [4] enrichplot_1.26.2      org.Mm.eg.db_3.20.0    clusterProfiler_4.14.3
##  [7] org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   IRanges_2.40.0        
## [10] S4Vectors_0.44.0       Biobase_2.66.0         BiocGenerics_0.52.0   
## [13] readr_2.1.5            dplyr_1.1.4            tidyr_1.3.1           
## [16] knitr_1.49            
## 
## loaded via a namespace (and not attached):
##  [1] DBI_1.2.3               gson_0.1.0              rlang_1.1.4            
##  [4] magrittr_2.0.3          DOSE_4.0.0              compiler_4.4.2         
##  [7] RSQLite_2.3.8           png_0.1-8               vctrs_0.6.5            
## [10] reshape2_1.4.4          stringr_1.5.1           pkgconfig_2.0.3        
## [13] crayon_1.5.3            fastmap_1.2.0           XVector_0.46.0         
## [16] labeling_0.4.3          utf8_1.2.4              rmarkdown_2.29         
## [19] tzdb_0.4.0              UCSC.utils_1.2.0        purrr_1.0.2            
## [22] bit_4.5.0               xfun_0.49               zlibbioc_1.52.0        
## [25] cachem_1.1.0            aplot_0.2.3             GenomeInfoDb_1.42.0    
## [28] jsonlite_1.8.9          blob_1.2.4              BiocParallel_1.40.0    
## [31] parallel_4.4.2          R6_2.5.1                bslib_0.8.0            
## [34] stringi_1.8.4           RColorBrewer_1.1-3      jquerylib_0.1.4        
## [37] GOSemSim_2.32.0         Rcpp_1.0.13-1           ggtangle_0.0.4         
## [40] R.utils_2.12.3          Matrix_1.7-1            splines_4.4.2          
## [43] igraph_2.1.1            tidyselect_1.2.1        qvalue_2.38.0          
## [46] yaml_2.3.10             codetools_0.2-20        lattice_0.22-6         
## [49] tibble_3.2.1            plyr_1.8.9              treeio_1.30.0          
## [52] withr_3.0.2             KEGGREST_1.46.0         evaluate_1.0.1         
## [55] gridGraphics_0.5-1      Biostrings_2.74.0       pillar_1.9.0           
## [58] ggtree_3.14.0           ggfun_0.1.7             generics_0.1.3         
## [61] vroom_1.6.5             hms_1.1.3               ggplot2_3.5.1          
## [64] munsell_0.5.1           scales_1.3.0            tidytree_0.4.6         
## [67] glue_1.8.0              lazyeval_0.2.2          tools_4.4.2            
## [70] data.table_1.16.2       babelgene_22.9          fs_1.6.5               
## [73] fastmatch_1.1-4         cowplot_1.1.3           grid_4.4.2             
## [76] ape_5.8                 colorspace_2.1-1        nlme_3.1-166           
## [79] GenomeInfoDbData_1.2.13 patchwork_1.3.0         cli_3.6.3              
## [82] fansi_1.0.6             gtable_0.3.6            R.methodsS3_1.8.2      
## [85] yulab.utils_0.1.8       sass_0.4.9              digest_0.6.37          
## [88] ggrepel_0.9.6           ggplotify_0.1.2         farver_2.1.2           
## [91] memoise_2.0.1           htmltools_0.5.8.1       R.oo_1.27.0            
## [94] lifecycle_1.0.4         httr_1.4.7              GO.db_3.20.0           
## [97] bit64_4.5.2
```
