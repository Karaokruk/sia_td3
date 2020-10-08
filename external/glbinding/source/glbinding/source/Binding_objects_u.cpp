
#include "Binding_pch.h"


using namespace gl;


namespace glbinding
{


Function<void, GLint, GLdouble> Binding::Uniform1d("glUniform1d");
Function<void, GLint, GLsizei, const GLdouble *> Binding::Uniform1dv("glUniform1dv");
Function<void, GLint, GLfloat> Binding::Uniform1f("glUniform1f");
Function<void, GLint, GLfloat> Binding::Uniform1fARB("glUniform1fARB");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform1fv("glUniform1fv");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform1fvARB("glUniform1fvARB");
Function<void, GLint, GLint> Binding::Uniform1i("glUniform1i");
Function<void, GLint, GLint64> Binding::Uniform1i64ARB("glUniform1i64ARB");
Function<void, GLint, GLint64EXT> Binding::Uniform1i64NV("glUniform1i64NV");
Function<void, GLint, GLsizei, const GLint64 *> Binding::Uniform1i64vARB("glUniform1i64vARB");
Function<void, GLint, GLsizei, const GLint64EXT *> Binding::Uniform1i64vNV("glUniform1i64vNV");
Function<void, GLint, GLint> Binding::Uniform1iARB("glUniform1iARB");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform1iv("glUniform1iv");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform1ivARB("glUniform1ivARB");
Function<void, GLint, GLuint> Binding::Uniform1ui("glUniform1ui");
Function<void, GLint, GLuint64> Binding::Uniform1ui64ARB("glUniform1ui64ARB");
Function<void, GLint, GLuint64EXT> Binding::Uniform1ui64NV("glUniform1ui64NV");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::Uniform1ui64vARB("glUniform1ui64vARB");
Function<void, GLint, GLsizei, const GLuint64EXT *> Binding::Uniform1ui64vNV("glUniform1ui64vNV");
Function<void, GLint, GLuint> Binding::Uniform1uiEXT("glUniform1uiEXT");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform1uiv("glUniform1uiv");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform1uivEXT("glUniform1uivEXT");
Function<void, GLint, GLdouble, GLdouble> Binding::Uniform2d("glUniform2d");
Function<void, GLint, GLsizei, const GLdouble *> Binding::Uniform2dv("glUniform2dv");
Function<void, GLint, GLfloat, GLfloat> Binding::Uniform2f("glUniform2f");
Function<void, GLint, GLfloat, GLfloat> Binding::Uniform2fARB("glUniform2fARB");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform2fv("glUniform2fv");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform2fvARB("glUniform2fvARB");
Function<void, GLint, GLint, GLint> Binding::Uniform2i("glUniform2i");
Function<void, GLint, GLint64, GLint64> Binding::Uniform2i64ARB("glUniform2i64ARB");
Function<void, GLint, GLint64EXT, GLint64EXT> Binding::Uniform2i64NV("glUniform2i64NV");
Function<void, GLint, GLsizei, const GLint64 *> Binding::Uniform2i64vARB("glUniform2i64vARB");
Function<void, GLint, GLsizei, const GLint64EXT *> Binding::Uniform2i64vNV("glUniform2i64vNV");
Function<void, GLint, GLint, GLint> Binding::Uniform2iARB("glUniform2iARB");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform2iv("glUniform2iv");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform2ivARB("glUniform2ivARB");
Function<void, GLint, GLuint, GLuint> Binding::Uniform2ui("glUniform2ui");
Function<void, GLint, GLuint64, GLuint64> Binding::Uniform2ui64ARB("glUniform2ui64ARB");
Function<void, GLint, GLuint64EXT, GLuint64EXT> Binding::Uniform2ui64NV("glUniform2ui64NV");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::Uniform2ui64vARB("glUniform2ui64vARB");
Function<void, GLint, GLsizei, const GLuint64EXT *> Binding::Uniform2ui64vNV("glUniform2ui64vNV");
Function<void, GLint, GLuint, GLuint> Binding::Uniform2uiEXT("glUniform2uiEXT");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform2uiv("glUniform2uiv");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform2uivEXT("glUniform2uivEXT");
Function<void, GLint, GLdouble, GLdouble, GLdouble> Binding::Uniform3d("glUniform3d");
Function<void, GLint, GLsizei, const GLdouble *> Binding::Uniform3dv("glUniform3dv");
Function<void, GLint, GLfloat, GLfloat, GLfloat> Binding::Uniform3f("glUniform3f");
Function<void, GLint, GLfloat, GLfloat, GLfloat> Binding::Uniform3fARB("glUniform3fARB");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform3fv("glUniform3fv");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform3fvARB("glUniform3fvARB");
Function<void, GLint, GLint, GLint, GLint> Binding::Uniform3i("glUniform3i");
Function<void, GLint, GLint64, GLint64, GLint64> Binding::Uniform3i64ARB("glUniform3i64ARB");
Function<void, GLint, GLint64EXT, GLint64EXT, GLint64EXT> Binding::Uniform3i64NV("glUniform3i64NV");
Function<void, GLint, GLsizei, const GLint64 *> Binding::Uniform3i64vARB("glUniform3i64vARB");
Function<void, GLint, GLsizei, const GLint64EXT *> Binding::Uniform3i64vNV("glUniform3i64vNV");
Function<void, GLint, GLint, GLint, GLint> Binding::Uniform3iARB("glUniform3iARB");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform3iv("glUniform3iv");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform3ivARB("glUniform3ivARB");
Function<void, GLint, GLuint, GLuint, GLuint> Binding::Uniform3ui("glUniform3ui");
Function<void, GLint, GLuint64, GLuint64, GLuint64> Binding::Uniform3ui64ARB("glUniform3ui64ARB");
Function<void, GLint, GLuint64EXT, GLuint64EXT, GLuint64EXT> Binding::Uniform3ui64NV("glUniform3ui64NV");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::Uniform3ui64vARB("glUniform3ui64vARB");
Function<void, GLint, GLsizei, const GLuint64EXT *> Binding::Uniform3ui64vNV("glUniform3ui64vNV");
Function<void, GLint, GLuint, GLuint, GLuint> Binding::Uniform3uiEXT("glUniform3uiEXT");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform3uiv("glUniform3uiv");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform3uivEXT("glUniform3uivEXT");
Function<void, GLint, GLdouble, GLdouble, GLdouble, GLdouble> Binding::Uniform4d("glUniform4d");
Function<void, GLint, GLsizei, const GLdouble *> Binding::Uniform4dv("glUniform4dv");
Function<void, GLint, GLfloat, GLfloat, GLfloat, GLfloat> Binding::Uniform4f("glUniform4f");
Function<void, GLint, GLfloat, GLfloat, GLfloat, GLfloat> Binding::Uniform4fARB("glUniform4fARB");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform4fv("glUniform4fv");
Function<void, GLint, GLsizei, const GLfloat *> Binding::Uniform4fvARB("glUniform4fvARB");
Function<void, GLint, GLint, GLint, GLint, GLint> Binding::Uniform4i("glUniform4i");
Function<void, GLint, GLint64, GLint64, GLint64, GLint64> Binding::Uniform4i64ARB("glUniform4i64ARB");
Function<void, GLint, GLint64EXT, GLint64EXT, GLint64EXT, GLint64EXT> Binding::Uniform4i64NV("glUniform4i64NV");
Function<void, GLint, GLsizei, const GLint64 *> Binding::Uniform4i64vARB("glUniform4i64vARB");
Function<void, GLint, GLsizei, const GLint64EXT *> Binding::Uniform4i64vNV("glUniform4i64vNV");
Function<void, GLint, GLint, GLint, GLint, GLint> Binding::Uniform4iARB("glUniform4iARB");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform4iv("glUniform4iv");
Function<void, GLint, GLsizei, const GLint *> Binding::Uniform4ivARB("glUniform4ivARB");
Function<void, GLint, GLuint, GLuint, GLuint, GLuint> Binding::Uniform4ui("glUniform4ui");
Function<void, GLint, GLuint64, GLuint64, GLuint64, GLuint64> Binding::Uniform4ui64ARB("glUniform4ui64ARB");
Function<void, GLint, GLuint64EXT, GLuint64EXT, GLuint64EXT, GLuint64EXT> Binding::Uniform4ui64NV("glUniform4ui64NV");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::Uniform4ui64vARB("glUniform4ui64vARB");
Function<void, GLint, GLsizei, const GLuint64EXT *> Binding::Uniform4ui64vNV("glUniform4ui64vNV");
Function<void, GLint, GLuint, GLuint, GLuint, GLuint> Binding::Uniform4uiEXT("glUniform4uiEXT");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform4uiv("glUniform4uiv");
Function<void, GLint, GLsizei, const GLuint *> Binding::Uniform4uivEXT("glUniform4uivEXT");
Function<void, GLuint, GLuint, GLuint> Binding::UniformBlockBinding("glUniformBlockBinding");
Function<void, GLuint, GLint, GLuint> Binding::UniformBufferEXT("glUniformBufferEXT");
Function<void, GLint, GLuint64> Binding::UniformHandleui64ARB("glUniformHandleui64ARB");
Function<void, GLint, GLuint64> Binding::UniformHandleui64NV("glUniformHandleui64NV");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::UniformHandleui64vARB("glUniformHandleui64vARB");
Function<void, GLint, GLsizei, const GLuint64 *> Binding::UniformHandleui64vNV("glUniformHandleui64vNV");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix2dv("glUniformMatrix2dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix2fv("glUniformMatrix2fv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix2fvARB("glUniformMatrix2fvARB");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix2x3dv("glUniformMatrix2x3dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix2x3fv("glUniformMatrix2x3fv");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix2x4dv("glUniformMatrix2x4dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix2x4fv("glUniformMatrix2x4fv");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix3dv("glUniformMatrix3dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix3fv("glUniformMatrix3fv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix3fvARB("glUniformMatrix3fvARB");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix3x2dv("glUniformMatrix3x2dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix3x2fv("glUniformMatrix3x2fv");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix3x4dv("glUniformMatrix3x4dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix3x4fv("glUniformMatrix3x4fv");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix4dv("glUniformMatrix4dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix4fv("glUniformMatrix4fv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix4fvARB("glUniformMatrix4fvARB");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix4x2dv("glUniformMatrix4x2dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix4x2fv("glUniformMatrix4x2fv");
Function<void, GLint, GLsizei, GLboolean, const GLdouble *> Binding::UniformMatrix4x3dv("glUniformMatrix4x3dv");
Function<void, GLint, GLsizei, GLboolean, const GLfloat *> Binding::UniformMatrix4x3fv("glUniformMatrix4x3fv");
Function<void, GLenum, GLsizei, const GLuint *> Binding::UniformSubroutinesuiv("glUniformSubroutinesuiv");
Function<void, GLint, GLuint64EXT> Binding::Uniformui64NV("glUniformui64NV");
Function<void, GLint, GLsizei, const GLuint64EXT *> Binding::Uniformui64vNV("glUniformui64vNV");
Function<void> Binding::UnlockArraysEXT("glUnlockArraysEXT");
Function<GLboolean, GLenum> Binding::UnmapBuffer("glUnmapBuffer");
Function<GLboolean, GLenum> Binding::UnmapBufferARB("glUnmapBufferARB");
Function<GLboolean, GLuint> Binding::UnmapNamedBuffer("glUnmapNamedBuffer");
Function<GLboolean, GLuint> Binding::UnmapNamedBufferEXT("glUnmapNamedBufferEXT");
Function<void, GLuint> Binding::UnmapObjectBufferATI("glUnmapObjectBufferATI");
Function<void, GLuint, GLint> Binding::UnmapTexture2DINTEL("glUnmapTexture2DINTEL");
Function<void, GLuint, GLuint, GLsizei, const void *, GLenum> Binding::UpdateObjectBufferATI("glUpdateObjectBufferATI");
Function<void, GLuint> Binding::UseProgram("glUseProgram");
Function<void, GLhandleARB> Binding::UseProgramObjectARB("glUseProgramObjectARB");
Function<void, GLuint, UseProgramStageMask, GLuint> Binding::UseProgramStages("glUseProgramStages");
Function<void, GLenum, GLuint> Binding::UseShaderProgramEXT("glUseShaderProgramEXT");



} // namespace glbinding
