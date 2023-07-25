pal_stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

pal_stallion_interp <- function(x = NULL) {
    if (is.null(x)) {
        return(colorRampPalette(pal_stallion))
    } else if (is.numeric(x) & length(x) == 1) {
        return(colorRampPalette(pal_stallion)(n)) 
    } else if (!is.null(x)) {
        return(colorRampPalette(pal_stallion)(length(unique(x))))
    }
}
