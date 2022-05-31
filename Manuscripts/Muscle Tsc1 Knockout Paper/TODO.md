# Pre-submission to-do list.

*1.* We need to complete a declaration of competing interests both in the manuscript file and on a separate form. 

*2.* The title page requires the full postal addresses of each affiliation and, if available, current email addresses for each author. 

*3.* On the title page, where your contact info is, they need a full postal address and telephone number. Please confirm the address I pulled from your UM facepage and add your office telephone number. 

*4.* It also specifies the affiliations be given as the one where the work was done, and new/current affiliations be included as a footnote. I don't know how strict on this they are, but since the writing/editing has happened at all our current workplaces, I'm sure we can leave the affiliations section as-is. 

*5.* The abstract guidelines are for a structured abstract with objective, methods, results and conclusions. I didn't change the abstract but was hoping we could get away with not changing it unless the journal explicitly asks. What are your thoughts? I'd hate to be triaged because of a formatting issue 

*6.* A highlights section is mandatory. I've thrown a few bullet points together as a separate file (a separate highlights file with 3-5 bullet points of max 85 characters per point). Please edit and add things as you see fit.

*7.* I added the keywords section. Please edit as you see fit.

*8.* An abbreviation section is required for any non-standard abreviations. I left it out, but add it if you think we need to. 

*9.* Manuscript sections should be numbered. i.e., 1. Introduction. And subsections should be numbered, i.e., 2. Methods, 2.1. Mice

*10.* If we want, we can include a graphical abstract.

*11.* Figures need to be in TIFF or JPEG with min 300dpi. Can also be PDF if vector drawings.

*12.* We have the option of converting our supplementary data into a data article, a new kind of article that houses and describes the data, known as Data in Brief. It would be given a DOI and be publically available upon publication. Mol Metab encourages researchers to submit data as a Data in Brief as an additional file upon submission of the manuscript revision.

#Things to keep in mind when preparing the submission cover letter
"Molecular Metabolism publishes hypothesis driven research generated with the highest standards that paves the way to a mechanistic understanding of energy homeostasis-associated behavior, physiology and dysfunction."






# First Submission ToDo List

Comments from MCM, DS, ARS, QT

No Comments from MP, KF

No Response from BL, JCH, JRR


## General Comments About Article

ARS: Very interesting study, but I think there are a lot of things missing that the reviewers will ding you on no matter where you send the paper.  I mentioned a few of these in the text.  I didn’t really do a serious rewrite of the paper since I think these missing experiments really need attention.  I think the biggest question is why there is such a big decrease in fat mass, with such a small effect on energy expenditure.  I get that there is a lot of muscle, but still…..  There are lots of missing things- same measurements in ND and HFD, incomplete reporting on metabolic cage studies, serum levels of insulin, glucose, ffas, leptin, cytokines, etc. It just needs a lot more.  I don’t think AJP is taking papers like this that aren’t really complete.  Also what about insulin sensitivity in muscle?  I think you also need measurements of fat and glucose oxidation- not that hard to do.  Sorry to sound so negative.  I think reviewers will be much harsher.

QT: I think Figure 4 and the text related to this figure need a bit of work.

## Comments to Think About

* **Both short term overfeeding and chronic obesity result in increased energy expenditure.**  - very controversial, I don’t agree.  Yes to short term, no to long term (ARS)
ES: As far as I know, the data support increased EE in obesity. Larger tissue mass=greater metabolic demand + decreased economy of movement/physical activity=greater TEE. Not sure what to think of this comment.

* **This is consistent with the hypothesis that the adiposity is inversely related to insulin sensitivity and that these mice are not lipodystrophic, but the mice with muscle Tsc1 knockout are protected from adipose tissue expansion via changes in energy balance.** Was this on normal or HFD?  Looks like ND.  How do you explain that?  You’d think that the muscle would be insulin resistant. (ARS).  
ES: Perhaps we need to include the chow clamp data, at least as a supplement, and then include something in the text to say there were no differences in insulin sensitivity for the ITT in chow fed mice and no differences in muscle glucose uptake during clamp, but increased glucose uptake in adipose under clamp conditions in chow mice. This would probably be enough. If there are fasted and refed mTSC1 muscle tissues/lysates somewhere, could always measure Akt phosphorylation or something we were asked to. 

Related to this maybe we should include the ITT and glucose clamp data showing no changes in insulin sensitivity on NCD (DB).
ES: Yes, and perhaps we should mention in the text that there was no difference in ITT/clamp on chow diet and  include the clamp data as a supplementary figure? Or we add it in as an early main figure. 

* For Mixed linear models.  I don’t think we should use Chi-square test here because they are continuous variables. I will check the R codes and let you know. (QT).
ES: I'm gonna leave that to QT to figure out.

* I think the biggest question is why there is such a big decrease in fat mass, with such a small effect on energy expenditure. (ARS).  
We could calculate out what the expected change in energy balance would be based on our measured intake/increases in energy expenditure and see if that comes close to matching.  It would be tricky because it wouldnt account for an adaptation, but for example if there was a 7% increase in EE how much less fat mass would be expected.  Ill think about this (DB)
ES: This is complicated. My thinking is that if this is easy enough to calculate, we should look at it. A 7 percent increase in EE every day, over the entire life of a mouse would cumulatively be quite a lot of excess energy expended. Otherwise, we could include it as a discussion point that we leave open for future investigation. That said, if we hand-wave about potential mechanisms (mTORC1-driven myokines, etc.), we could open ourselves up for a ton of additional experiment requests. We have the FGF21 data, but I don't know if that would be enough.

* incomplete reporting on metabolic cage studies, serum levels of insulin, glucose, ffas, leptin, cytokines, etc. It just needs a lot more. (ARS).  
Not sure what is missing about metabolic cage stuff
ES: Maybe it's because we are only showing heat values (which is correct for TEE), as opposed to also including VO2 and VCO2 data (from which heat is calculated). Some people like to see everything so you can get a sense of whole-body coupling efficiency. Perhaps we should include the regression curves for VO2/VCO2 and lean mass. I personally like to see this data presented in papers.

Im not sure if we have insulin/FFAs/leptin/cytokines but i can look into what 70d old blood we have available.  
ES: We have insulin (there is a figure in my copy of the clamp data folder) and fasting glucose (itt folder) for chow. We also have serum FGF21. No NEFA or TG, as far as I can remember (there's no evidence this was done, but liver TG's was measured). Maybe we could put a table in the paper with what we do have? We didn't do a bioplex assay on the mTSC1 serum, but if there is serum left it might be something to consider doing. 

* I think you also need measurements of fat and glucose oxidation- not that hard to do (ARS).  
ES: As we already discussed, we can do this on a whole-body level using the CLAMS data. Measuring these things in muscle is not difficult but would require new mice. Kinda related but not really, we also have muscle glycogen (increased in fasted mice compared to CON's) and muscle TG (I can only find fasted data, but there was no difference). 


## Things to Change

* Can we stick with the same color scheme used in Figure 1, i.e. knockout is in cyan (greenish color) and wild-type/control is in pink/orange. We have these 2 colors switched in Figure 2. (QT).
ES: Agree. I think the inconsistencies in color scheme between figures is something I mentioned before we sent the paper out. Might be difficult where we have males and females on the same graph.
* Same here. Can we use color for the bar plot just like the one in Figure 1? It is easier if we use the same color scheme throughout so when we look at the figures, we know cyan is knockout and red/pink/oranges (pick one) is control.  If not, then everything should be in dark gray or light gray with labels (like in the figure 3C). Too many colors are used here. (QT)
ES: Agree.
* Add Body weight to Figure 3.
ES: Ok. Sure, this could go here. Otherwise we could report end-point weights in the text.
* Add insulin sensitivity to Figure 2, including an explanation of why they are **not** insulin resistant.
 ES: Actually, this would probably be a good place to include the clamp data. Or as a separate figure between current Fig's 2 and 3. 
* Rewrite section describing Figure 4 or the figure, so that panels are described in the same order that they appear (see QT comments).
ES: I have re-arranged some of this stuff, including in the figure (made and alternate figure 4 in case it doesn't work for others).

## For Resubmission

* Add relationship between TDEE and Ageing