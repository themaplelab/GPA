void swap(char **p, char **q);

int main(){
    char A[1], B[1];
    char *a, *b, **t1, **t2;
    t1 = &a;
    t2 = &b;
    a = &A;
    b = &B;

    if(A[0]){
        *a = 'A';
    }
    else{
        swap(t1, t2);
        *a = 'B';
    }

    *a = '?';
}

void swap(char **p, char **q){
    /*
    
    %0 = load i8*, i8** %p, align 8, !tbaa !13
    %1 = load i8*, i8** %q, align 8, !tbaa !13
    store i8* %1, i8** %p, align 8, !tbaa !13
    store i8* %0, i8** %q, align 8, !tbaa !13
    
    */
    char *tmp = *p;
    char **x = p;
    char *a = *q;
    *p = a;
    *q = tmp;

}