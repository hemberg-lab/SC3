# See https://github.com/ropensci/RSelenium/issues/67 
# for explanation of why this function exists.
stopSeleniumServer <- function()
{
    if (.Platform$OS.type == "windows")
    {
        procs <- system2("wmic",
            "path win32_process get Caption,Processid,Commandline",
            stdout=TRUE, stderr=NULL)
        idx <- grep("selenium-server-standalone.jar", procs)
        if (!length(idx)) # selenium not running?
            return()
        proc <- procs[idx]
        proc <- trimws(proc)
        segs <- strsplit(proc, " ", TRUE)[[1]]
        pid <- segs[length(segs)]
        system2("taskkill", c("/pid", pid, "/f"), stdout=NULL, stderr=NULL)
    } else {
        # Unfortunately field width specifiers do not work on Mac
        # so we need two commands to get pid and command non-truncated.
        pidcmd <- system2("ps", "wwaxo pid,command", stdout=TRUE, stderr=NULL)
        pidusr <- system2("ps", "wwaxo pid,user", stdout=TRUE, stderr=NULL)
        pidusr <- trimws(pidusr)
        whoami <- system2("whoami", stdout=TRUE, stderr=NULL)
        pididx <- grep(paste0(whoami, "$"), pidusr)
        pidusr <- pidusr[pididx]
        tmp <- strsplit(pidusr, " ")
        pids <- unlist(lapply(tmp, function(x) x[1]))
        procs_idx <- grep("selenium-server-standalone.jar", pidcmd)
        for (proc in pidcmd[procs_idx])
        {
            pid <- strsplit(proc, " ")[[1]][1]
            if (pid %in% pids) 
            {
                system2("kill", c("-9", pid))
            }
        }
    }
}