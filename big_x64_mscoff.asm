;	Copyright (c) 2012 Christopher A. Taylor.  All rights reserved.
;
;	Redistribution and use in source and binary forms, with or without
;	modification, are permitted provided that the following conditions are met:
;
;	* Redistributions of source code must retain the above copyright notice,
;	  this list of conditions and the following disclaimer.
;	* Redistributions in binary form must reproduce the above copyright notice,
;	  this list of conditions and the following disclaimer in the documentation
;	  and/or other materials provided with the distribution.
;	* Neither the name of LibCat nor the names of its contributors may be used
;	  to endorse or promote products derived from this software without
;	  specific prior written permission.
;
;	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
;	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
;	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;	ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
;	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
;	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
;	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
;	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
;	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;	POSSIBILITY OF SUCH DAMAGE.

; Assembly code written for FASM (flatassembler.net)

format MS64 COFF

section '.code' code readable executable 


;-----------------------------------------------------------------------------
;                           Function Declarations
;-----------------------------------------------------------------------------
public _mul_mod64 as '_mul_mod64'

;-----------------------------------------------------------------------------
;                                 mulmod64
;-----------------------------------------------------------------------------
; Arguments:
;   rcx      (arg 1) : A
;   rdx      (arg 2) : B
;   r8       (arg 3) : Modulus
; Output:
;   return A * B % M
;-----------------------------------------------------------------------------
_mul_mod64:

	mov rax, rcx
	mul rdx
	div	r8
    mov rax, rdx
    ret


;-----------------------------------------------------------------------------
;                    Appendix: x86-64 Calling Convention
;-----------------------------------------------------------------------------
; Callee save: rsp(stack), rbx(base), rbp(frame), r12-r15, rdi, rsi
; Caller save: rax(ret), rcx, rdx, r8-r11
