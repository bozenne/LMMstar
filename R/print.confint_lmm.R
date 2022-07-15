### print.confint_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 15 2022 (17:48) 
## Version: 
## Last-Updated: jul 15 2022 (18:02) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.confint_lmm
##' @export
print.confint_lmm <- function(x, digit = 3, detail = FALSE, ...){
    print(as.data.frame(x), digits = digit)
    message.backtransform <- attr(x,"back-transform")

    if(detail>0){

        cat("\n   Note: ")

        method.multiplicity <- attr(x,"method")
        add.space <- ""

        if(!is.null(method.multiplicity) && method.multiplicity!="none"){
            
            if(NROW(x)==1){
                short2text <- stats::setNames(c("confidence interval","confidence interval","p-value"),c("lower","upper","p.value"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }else{
                short2text <- stats::setNames(c("confidence intervals","confidence intervals","p-values"),c("lower","upper","p.value"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }
            if(method.multiplicity == "bonferroni"){
                txt.method <- "Bonferroni"
            }else if(method.multiplicity %in% c("single-step","single-step2")){
                txt.method <- "max-test adjustment"
            }else{
                txt.method <- method.multiplicity
            }
            add.space <- "         "
            cat(paste(txt,collapse = ", ")," have been adjusted for multiplicity using ",method.multiplicity,". \n",sep="")
        }
        
        if(!is.null(message.backtransform) && any(!is.na(message.backtransform$FUN))){
            message.backtransform <- message.backtransform[!is.na(message.backtransform$FUN),,drop=FALSE]

            if(any(message.backtransform[,setdiff(names(message.backtransform), "FUN")] == FALSE)){
                warning("Could not back-transform everything.\n")
            }

            if(NROW(x)==1){
                short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }else{
                short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }
            
            cat(add.space,paste(txt,collapse = ", ")," have been back-transformed",sep="")
            if(detail>=0.5){
                cat(" (",paste0(paste(rownames(message.backtransform),collapse = "/")," parameters with ",paste(message.backtransform$FUN,collapse="/")),"). \n", sep ="")
            }
            cat(".\n")
        }
    }
    return(invisible(NULL))
}


##----------------------------------------------------------------------
### print.confint_lmm.R ends here
