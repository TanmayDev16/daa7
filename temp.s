nasm
; Pair class template in Assembly

section .data
    pair_open_paren db '(', 0
    pair_comma db ',', 0
    pair_close_paren db ')', 0

section .text
    global Pair_constructor
    global Pair_make_pair
    global Pair_print
    global Pair_get_first
    global Pair_get_second

; Default constructor
Pair_constructor:
    push ebp
    mov ebp, esp
    ; Call default constructors for first and second members
    ; ...
    pop ebp
    ret

; Parameterized constructor
Pair_constructor_params:
    push ebp
    mov ebp, esp
    ; Store first and second parameters in object
    ; ...
    pop ebp
    ret

; make_pair function
Pair_make_pair:
    push ebp
    mov ebp, esp
    ; Store first and second parameters in object
    ; ...
    pop ebp
    ret

; print function
Pair_print:
    push ebp
    mov ebp, esp
    ; Print open parenthesis
    push pair_open_paren
    call print_string
    add esp, 4
    ; Call print function for first member
    ; ...
    ; Print comma
    push pair_comma
    call print_string
    add esp, 4
    ; Call print function for second member
    ; ...
    ; Print close parenthesis
    push pair_close_paren
    call print_string
    add esp, 4
    pop ebp
    ret

; get_first function
Pair_get_first:
    push ebp
    mov ebp, esp
    ; Return first member
    ; ...
    pop ebp
    ret

; get_second function
Pair_get_second:
    push ebp
    mov ebp, esp
    ; Return second member
    ; ...
    pop ebp
    ret

; Helper function to print a string
print_string:
    push ebp
    mov ebp, esp
    mov eax, [ebp+8]  ; Get the string pointer
    mov ebx, eax      ; Save the string pointer
print_loop:
    mov al, [eax]     ; Load the next character
    cmp al, 0         ; Check if it's the null terminator
    je print_done     ; If so, exit the loop
    push eax          ; Push the character on the stack
    call print_char   ; Call print_char to print it
    pop eax           ; Pop the character off the stack
    inc eax           ; Move to the next character
    jmp print_loop    ; Loop back
print_done:
    pop ebp
    ret

print_char:
    ; Implementation to print a single character
    ; ...
    ret

