# Needed Packages:
package(limma) # from vienna package

# normal call of script - THIS IS MANDATORY TO DO IT AT LEAST ONCE SINCE THE SCRIPT IS CREATING DATA-FILES!
analyze()

# subplot of a defined region
analyze(whole=FALSE,s=[start coordinate],e=[end coordinate])
